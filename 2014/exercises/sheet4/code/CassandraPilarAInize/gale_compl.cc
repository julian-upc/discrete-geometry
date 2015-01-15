/* Copyright (c) 2014
 Cassandra Pilar and Ainize FME/UPC 2014
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2, or (at your option) any
 later version: http://www.gnu.org/licenses/gpl.txt.
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 --------------------------------------------------------------------------------
 */

#include "polymake/Set.h"
#include "polymake/Vector.h"
#include "polymake/SparseVector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/SparseMatrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/PowerSet.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"
#include "polymake/Integer.h"
#include "polymake/Map.h"
#include "polymake/client.h"
#include <string>
#include <iostream>
#include <set>
#include <vector>
#include <list>
#include <iterator>


namespace polymake { namespace polytope {
    
    Matrix<Rational> enumerate_configurations(int e,int n,int m)
    {
        Matrix<Rational> equations(e+1,e*n+1);//we needd to add the equatin v00=m at the end
        for (int j=0;j<e;++j){
            for (int i=0;i<n;++i){
                equations(j,i*e+j+1)=1;}}
        equations(e,0)=-m;
        equations(e,1)=1;
        
        int lign_counter = 0;
        Matrix<Rational> ineq(e+2*e*n+(n*(n-1))/2,e*n+1);
        for(int i=0;i<e-1;++i){ //leave out last vector
            ineq[lign_counter][i+1]=1;
            ineq[lign_counter][i+2]=-1;
            ++lign_counter;
        }
        ineq[lign_counter][e]=1;
        ++lign_counter;
        
        for(int i=0;i<n-1;++i) {
            for(int j=i+1;j<n;++j) {
                ineq[lign_counter][e*i+1]=1;
                ineq[lign_counter][e*j+1]=-1;
                ++lign_counter;
            }
        }
        
        
        for (int i=1;i<e*n+1;++i){
            ineq[lign_counter][0]=m;
            ineq[lign_counter][i]=1;
            ++lign_counter;
        }
        for (int i=1;i<e*n+1;++i){
            ineq[lign_counter][0]=m;
            ineq[lign_counter][i]=-1;
            ++lign_counter;
        }
        
        perl::Object P("Polytope");
        P.take("EQUATIONS") << equations;
        P.take("INEQUALITIES")<< ineq;
        const Matrix<Rational> C = P.CallPolymakeMethod("LATTICE_POINTS");
        
        return C;
        
        
    }
    
    
    class Galepolytope{
    public:
        Vector<Integer> vec;
        Matrix<Rational> G;
        Matrix<Integer> C;//circuits
        int n_facets;
        Rational maxn_vertex;
    } ;
    
    
    
    bool lexicoordered(const Vector<Integer>& V ,int n, int e){
        for (int i=1; i<n*e-1; i += e) {
            if (V.slice(i,e) < V.slice(e+i, e)) {
                return false;
            }
        }
        return true;
    }
    
    bool cocircuits_or_fail(const Matrix<Rational>& G,
                            std::vector<Set<int> >& positive_parts,
                            std::vector<Set<int> >& negative_parts,
                            bool allow_fail=true)
    {
        //Enumerate all the cocircuits of G
        const int n = G.rows(), e = G.cols();
        for (Entire<Subsets_of_k<const sequence&> >::const_iterator sit = entire(all_subsets_of_k(sequence(0, n), e-1)); !sit.at_end(); ++sit) {
            const Matrix<Rational> ker = null_space(G.minor(*sit, All));
            if (ker.rows() >= 2){
                continue;}
            //Evaluate all the points in the ker to see if they are positive or negative (0 if they are in the defining subset)
            Set<int> plus, minus;
            for (int i=0; i<n; ++i) {
                const Rational val(ker[0] * G[i]);
                if (val > 0) plus += i;
                else if (val < 0) minus += i;
            }
            if (allow_fail &&
                (plus.size() == 1 ||
                 minus.size() == 1 ||
                 plus.size() == 0 && minus.size() == 0)) // not full-dimensional
                return false;
            positive_parts.push_back(plus);
            negative_parts.push_back(minus);
        }
        return
        positive_parts.size() > 0 ||//check how many potive and negative points we get
        negative_parts.size() > 0;
    }
    
    template<typename E>
    bool compare_matrices(const Matrix<E>& M1,const Matrix<E>& M2)
    {
        int nequal=0;
        for (int i=0;i<M1.rows();++i){
            for (int j=0; j<M2.rows();++j){
                if (M1[i]==M2[j]) {
                    nequal++;
                }
            }
        }
        return nequal == M1.rows();
    }
    
    bool compare_polytopes(const Galepolytope& P1, const Galepolytope& P2)
    {
        if (P1.n_facets<P2.n_facets) {
            return false;
        } else if (P1.n_facets == P2.n_facets) {
            if (P1.maxn_vertex < P2.maxn_vertex) {
                return false;
            }
        }
        return true;
    }
    
    //empty list of polytopes
    std::list<Galepolytope> polylist;
    std::list<Matrix<Rational> > galelist;
    
    std::list<Matrix<Rational> > gale_compl(int e,int n,int m)
    {
        const Matrix<Rational> configs=enumerate_configurations(e,n,m);
        //for each configuration:
        
        for(Entire <Rows< Matrix<Rational> > >::const_iterator rit = entire(rows(configs)); !rit.at_end(); ++rit){
            const Vector<Integer> v(*rit);
            // first check if it is ordered
            if (!lexicoordered(v,n,e)) {
                continue;
            }
            
            
            Matrix<Rational> G(n, e);// you need to initialize the matrix with its dimensions
            for (int j=1;j<e+1;++j){
                G[0][j-1]=v[j];}
            for (int i=1; i<n; ++i) {
                for (int j=0; j<e; ++j) {
                    G[i][j] = v[i*e+j];
                }
            }
            
            //check if convex with cocitcuits_or_fail
            std::vector<Set<int> > positive_parts;
            std::vector<Set<int> > negative_parts;
            if (!cocircuits_or_fail(G, positive_parts, negative_parts)) {
                continue;
                
            }
    
            
            //circuits
            //Matrix<Rational>  ker = ones_matrix<Rational>(1, n);
            Matrix<Rational>  ker=null_space(T(G));
            
            std::cout << "rows of the kernel  "<< ker.rows() <<endl;
            
            
            Set<int>  positive_index;
            for (int i=0; i<ker.rows(); i++) positive_index.insert(i);
            //Take the sign
            Matrix<Rational> circuits=ker;
            for (int i=0; i<ker.rows();++i) {
                for (int j=0; j<ker.cols();++j) {
                    if (ker[i][j]>0) {
                        circuits[i][j]=1;
                    } else if (ker[i][j]<0) {
                        circuits[i][j]=-1;
                        positive_index.erase(i);
                    }
                }
            }
            std::cout<< "pos index size  "<< positive_index.size()<<endl;
            circuits=circuits.minor(positive_index,All);
            
            Rational max=0;
            for (int i=0; i<circuits.rows(); ++i) {
                Rational sum=0;
                for (int j=0; j<circuits[i].size(); ++j) {
                    sum += circuits[i][j];
                }
                if (sum>max){
                    max=sum;
                }
            }
            
            
            //P.maxn_vertex=max;
            //Galepolytope Po(v,G,circuits,circuits.rows(),max);
            Galepolytope Po;
            Po.vec=v;
            Po.G=G;
            Po.C=circuits;
            Po.n_facets=circuits.rows();
            Po.maxn_vertex=max;
            
            
            //add to list
            
            polylist.push_back(Po);
        }
        
        //analyze the list we got:
        //order the list by number of facets:
        std::cout<< "size of list"<< polylist.size()<<endl;
        
        polylist.sort(compare_polytopes);
        
        int t=0;
        for (std::list<Galepolytope>::const_iterator iterator = polylist.begin(), end = polylist.end(); ++iterator != end; ++iterator)
        {   if (++iterator==end){
            Galepolytope P1=*iterator;
            Galepolytope P2=*++iterator;
            if (P1.n_facets != P2.n_facets)
            {galelist.push_back(P1.G);
                galelist.push_back(P2.G);
                t=t+2;
                continue;}
            if (P1.maxn_vertex != P2.maxn_vertex)
            {galelist.push_back(P1.G);
                galelist.push_back(P2.G);
                t=t+2;
                continue;}
            if( compare_matrices(P1.C,P2.C)) {
                galelist.push_back(P1.G);
                t=t+1;
                continue;}
        }
        else {
            Galepolytope P1=*iterator;
            Galepolytope P2=*++iterator;
            if (P1.n_facets != P2.n_facets)
            {galelist.push_back(P1.G);
                t=t+1;
                continue;}
            if (P1.maxn_vertex != P2.maxn_vertex)
            {galelist.push_back(P1.G);
                t=t+1;
                continue;}
            if( compare_matrices(P1.C,P2.C)) {continue;}
        }}
        return galelist;
        
    }
    UserFunctionTemplate4perl("# @category Computations"
                              "# Enumerate , up to combinarotial equivalence, all balanced configurations with integral Gale complexity equal to m."
                              "# @param Int e Dimension of the gale diagram vector space"
                              "# @param Int n Number of vertices of the Gale diagram"
                              "# @param Int m is the integral Gale complexity"
                              "# @return list<Matrix<Rational>> The list of all the matrices corresponding to the Gale diagrams.",
                              "gale_compl($$$)" );
    
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End: