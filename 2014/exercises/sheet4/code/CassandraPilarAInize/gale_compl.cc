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



namespace polymake { namespace polytope {

Matrix<Rational> enumerate_configurations(int e,int n,int m)
{
    Matrix<Rational> equations(e+1,e*n+1);//we needd to add the equatin v00=m at the end
    for (int j=0;j<e;++j){
        for (int i=0;i<n;++i){
            equations(j,i*e+j+1)=1;}}
    equations(e,0)=-m;
    equations(e,1)=1;
    
    ListMatrix<SparseVector<Rational > > inequalities(0,e*n+1);
    for(int i=0;i<e-1;++i){ //leave out last vector
        SparseVector<Rational>ineq(e*n+1);
        ineq[i+1]=1;
        ineq[i+2]=-1;
        inequalities/=ineq;
    }
    SparseVector<Rational>ineq(e*n+1);
    ineq[e]=1;
    inequalities/=ineq;
    for(int i=0;i<n-1;++i) {
        for(int j=0;j<e;++j) {
            SparseVector<Rational>ineq(n*e+1);
            ineq[e*i+j+1]=1;
            ineq[(e+1)*i*j+1]=-1;
            inequalities /= ineq;
        }
    }
    
    perl::Object P("Polytope");
    P.take("EQUALITIES");
    P.take("INEQUALITIES");
    const Matrix<Rational> C = P.CallPolymakeMethod("LATTICE_POINTS");
    
    return C;
}

    class Configuration{
        Vector<Integer> vec;
        Matrix<Rational> G;
        Matrix<Integer> C;//cocircuits
        
        Configuration(const Vector<Integer>& _v,
                      const Matrix<Rational>& _G,
                      const Matrix<Integer>& _C)
        : vec(_v)
        , G(_G)
        , C(_C) {}
    };
    
    // Configuration conf(vec1, G1, C1);
    // std::vector<Configuration> conf_vec;
    // conf_vec.push_back(conf);
    // w += conf.vec:

Matrix <Rational> gale_compl(int e,int n,int m)
{
    const Matrix<Rational> configs=enumerate_configurations(e,n,m);
    //interating through all the conf we got for validation:
    //first chech if it is ordered
    for(Entire <Rows<Matrix<Integer> > >::const_iterator rit=entire(rows(configs)); !rit.at_end(); ++rit){
        const Vector<Integer> configs[*rit];
         for (int i=0; i<n; i += e+1) {
            if ( lexicographical_compare(configs[*rit][0]+i,
                                              configs[*rit][0]+e+i,
                                              configs[*rit][0]+e+1+i,
                                              configs[*rit][0]+2*e+1+i)==0){
                ignore(configs[*rit]);}}// ignore the not ordered configs. (the iterator rit is not correctly defined)
             //Now enumerate cocircuits:
             
        
    }
}
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    


UserFunctionTemplate4perl("# @category Computations"
                          "# Compute a vector that selects the inclusion-minimal face "
			  "# that contains the vertices indexed by I."
                          "# @param Polytope P the input polytope"
                          "# @param Set I the indices of some vertices of P"
                          "# @return Vector",
			  "face_selector(Polytope, Set)" );

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
