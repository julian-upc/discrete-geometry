//Darrera edició: 12 de gener de 2015

#include <iostream>
#include <vector>
#include <cmath>
#include <set>
#include "polymake/Vector.h"
#include "polymake/SparseVector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Array.h"
#include "polymake/linalg.h"
#include <math.h>
#include <set> 
#include <algorithm> 

using namespace std;

namespace polymake{ namespace std{

typedef Matrix <Rational> MR;
// typedef ListMatrix <SparseVector<Rational> > LMR;
typedef vector< vector< vector<int> > > Tri_vector;
typedef vector< vector<int> >  Matriz;



struct gale_conf{
	Matrix<Rational> vectors; //to be able to calculate kernel breaks the configuration into vectors

	Matrix <Rational> transform; //stores the transformation of the gale configuration

	set<vector<int> > cocircuits; //cocircuits of the transform

	set<vector<int> > facets; //stores the facet lattice cocircuits of the transform
};


bool are_different(  vector<int> v1,  vector<int> v2) { 
	int minim;
	minim =  min(v1.size(), v2.size());
	for (int i = 0; i < minim; ++i) {
		if (v1[i] != v2[i]) return true;
	}
	return false;
}


int compare( int a, int b, int c ) {
    if ( a <= 0 || b <= 0 || c < 0 ) {
        throw invalid_argument( "Received invalid values" );
    }
    return false;
}


// devuelve matriz con k columnas con el indice en cada
vector< vector<int> > combinations(int n, int k) {
  
   vector<int> combi_actual;
    vector<  vector<int> > combis;
    vector<bool> v(n);
    fill(v.begin() + k, v.end(), true);
 
   do {
       for (int i = 0; i < n; ++i) {
           if (!v[i]) {
               combi_actual.push_back(i);
           }
       }
       combis.push_back(combi_actual);
       combi_actual.clear();
   } while ( next_permutation(v.begin(), v.end()));
   return combis;
}


set < vector<int> > generate_cocircuits ( const Matrix <Rational> & M, const  int n, const  int e){
	Matriz c = combinations (n, e);
	set<vector<int> > cocircuits;
	for (unsigned int i=0; i<c.size(); ++i){
		Vector <Rational>  x (e);
		Vector <Rational>  b (e,0);
		MR A (e,e);
		for (int j=0; j<n; ++j){
			A[j]=M[c[i][j]];
		}
		x=lin_solve (A,b);
		//find the cocircuit defined by this hiperplane
		Vector <Rational> y =x*M;
		vector < int> z (n);
		for (int j=0; j<n; ++j){
			if (y[j]>0) {
				z[j]=+1;
			}
			else if (y[j]==0){
				z[j]=0;
			} 
			else if (y[j]<0){
				z[j]=-1;
			} 
		}
		cocircuits.insert(z);
	}
	return cocircuits;
}


set<vector< int> > generate_cocircuits_trans ( const Matrix <Rational> & M, const  int n, const  int d){//d=n-e-1
	Matriz c = combinations (n, d);
	set<vector < int> > cocircuits;
	for (unsigned int i=0; i<c.size(); ++i){
		Vector <Rational>  x (d);
		Vector <Rational>  b (d,0);
		MR A (d,d);
		for (int j=0; j<n; ++j){
			for(int k=0; k<d; ++k){
				A[j][k]=M[c[i][j]][k+1];
			}
		}
		x=lin_solve (A,b);
		//find the cocircuit defined by this hiperplane
		Vector <Rational> x2 (d+1,1);
		for (int j=0; j<d; ++j){
			x2[j+1]=x[j];
		}
		Vector <Rational> y =x2*M;
		vector < int> z (n);
		for (int j=0; j<n; ++j){
			if (y[j]>0) {
				z[j]=+1;
			}
			else if (y[j]==0){
				z[j]=0;
			} 
			else if (y[j]<0){
				z[j]=-1;
			} 
		}
		cocircuits.insert(z);
	}
	return cocircuits;
}



bool compare_cocircuits ( set<vector< int> > & cocircuits1, set<vector< int> > & cocircuits2){
	if (cocircuits1!=cocircuits2) return false;
	return true;
}



set <vector <int> > facet_lattice (const set<vector< int> > & cocircuits){
	set <vector <int> > facets;
	for (set<vector< int> >::iterator it= cocircuits.begin(); it !=cocircuits.end(); ++it){
		int sign_facet=0;
		bool is_facet=true;
		vector<int> v=*it;
		for (unsigned int i=0; i<v.size() and is_facet; ++i){
			if (v[i]!=0 and sign_facet==0){
				sign_facet=v[i];
			}
			if (v[i]>0 and sign_facet<0){
				is_facet=false;
			}
			if (v[i]<0 and sign_facet>0){
				is_facet=false;
			}
		}
		if (is_facet){
			for (unsigned int i=0; i<v.size() and is_facet; ++i){//with this all the facets are positive
				v[i]=sign_facet*v[i];
			}
			facets.insert(v);
		}
	}
	return facets;
}




vector <gale_conf> galecomplexity (const  int e, const   int n, const  int m) {
 	int kmax, gcompl, val_actual, val_0;// no la usas! Gsize;
    	set<  vector<int> > Coor1D;
	vector<int> coor_actual;
	vector<bool> lexi_actual;
	vector<  vector<int> > conf_actual, conf_0;
	vector<  vector<bool> > Gale_lexi, Gale_lexi0; 
	Tri_vector Sk_1, Sk, Gale_conf, Gale_0;
	bool is_lexi;
    	compare(e, n, m);
    	kmax = floor(n/2)*m;
    	// Generación de configuraciones unidimensionales
    	// inicializacion
    	if (n > 1) {
    		Sk_1.resize(n-1);
    		Sk.resize(n-1);
    		for (int g = 0; g < n-1; ++g) {
      			Sk_1[g].resize(1);
      			Sk_1[g][0].resize(g+1);
    		}
    	// emparejado primero
    		for (int g = 0; g < 1; ++g) {
      			gcompl = n - 2 - g;
      			coor_actual = Sk_1[g][0];
      			coor_actual.insert(coor_actual.end(),Sk_1[gcompl][0].begin(),Sk_1[gcompl][0].end()); 
      			Coor1D.insert(coor_actual);
    		}
  	}
  	else {
    		coor_actual.clear();
    		coor_actual.resize(1);
    	C	oor1D.insert(coor_actual);
  	}
  	//write_mat(Coor1D);


  
  	for (int k = 1; k < kmax + 1; ++k) {
    		for (int g = 0; g < n-1; ++g) { // el nuemro de puntos en la config
      			for (unsigned int i = 0; i < Sk_1[g].size(); ++i) { // el numero de opciones 
				val_0 = Sk_1[g][i][0];	
				for (int j = 0; j < g + 1; ++j) { // las coordenadas, la que aumentar
	  				val_actual = Sk_1[g][i][j];
	  				if ( (j == 0) && (val_actual < m) ) {
	    				coor_actual = Sk_1[g][i];
	    				coor_actual[0] += 1;
	    				Sk[g].push_back(coor_actual);
          			}
	  				else if ( ((val_0 - val_actual) == 1) && (val_actual < m) && ((Sk_1[g][i][j-1] - val_actual) == 1) ){
	    				coor_actual = Sk_1[g][i];
	    				coor_actual[j] += 1;
	    				Sk[g].push_back(coor_actual);
          			}
				}
      			}
    		}
   		// emparejado 
    		for (int g = 0; g < ceil((n-1)/2); ++g) {
      			gcompl = n - 2 - g;
      			for (unsigned int i = 0; i < Sk[g].size(); ++i) {
        			for (unsigned int j = 0; j < Sk[gcompl].size(); ++j) {
  	  				coor_actual = Sk[g][i];
	  				for (unsigned int l = 0; l < coor_actual.size(); ++l) {
	    				coor_actual[l] *= -1;
	  				}
	  				coor_actual.insert(coor_actual.end(),Sk[gcompl][j].begin(),Sk[gcompl][j].end());
	   				sort(coor_actual.rbegin(), coor_actual.rend()); // reverse iterators
	  				//Coor1D.push_back(coor_actual);
					Coor1D.insert(coor_actual);	  
	  				coor_actual = Sk[gcompl][j];
	  				if (are_different(coor_actual, Sk[g][i])) {
	    					for (unsigned int l = 0; l < coor_actual.size(); ++l) {
	      						coor_actual[l] *= -1;
	    					}
	    					coor_actual.insert(coor_actual.end(),Sk[g][i].begin(),Sk[g][i].end());
	    					sort(coor_actual.rbegin(), coor_actual.rend()); // reverse iterators
	    					Coor1D.insert(coor_actual);
	  				}
        			}
      			}
    		}
    		Sk_1 = Sk; 
    		Sk.clear();
    		Sk.resize(n-1);
  	}
  
  	// Primera dimension
  	conf_actual.resize(n);
  	lexi_actual.resize(n);
  	for (int i = 0; i < n; ++i) {
    		conf_actual[i].resize(e);
  	}
  	for ( set<  vector<int> >::iterator it = Coor1D.begin(); it != Coor1D.end(); ++it) {
    		coor_actual = (*it);
    		if (coor_actual[0] == m) {
      			for (int i = 0; i < n; ++i) {
   					conf_actual[i][0] = coor_actual[i];
					if (i == 0) {
	  					lexi_actual[0] = true; 
					}      
					else if ((coor_actual[i-1] - coor_actual[i]) > 0) {
	  					lexi_actual[i] = true;
					}
					else {
		  				lexi_actual[i] = false;
					}
      			}
      			Gale_conf.push_back(conf_actual);
      			Gale_lexi.push_back(lexi_actual);
    		}   
    	//conf_actual.clear();
  	}
  	// Ampliación hasta e dimensiones
  	for (int d = 1; d < e; ++d) {
    		Gale_0 = Gale_conf;
    		Gale_lexi0 = Gale_lexi;
    		Gale_conf.clear();
    		Gale_lexi.clear();
    		for ( set<  vector<int> >::iterator it = Coor1D.begin(); it != Coor1D.end(); ++it) {
      			coor_actual = (*it);
      			for (unsigned int i = 0; i < Gale_0.size(); ++i) {
        			conf_0 = Gale_0[i];
        			lexi_actual = Gale_lexi0[i];
        			do {
 	  				is_lexi = true;
          				if (conf_0[0][d-1] < coor_actual[0]) {
	    					is_lexi = false;
	    					continue;
	  				}
	  				for (int j = 1; j < n; ++j) {
            					if ( (!lexi_actual[j]) && ((coor_actual[j-1] - coor_actual[j]) < 0) ) {
	      						is_lexi = false;
	      						break;
	    					}
	    					else if ( (!lexi_actual[j]) && ((coor_actual[j-1] - coor_actual[j]) > 0) ) {
	      						lexi_actual[j] = true;
  	    					}
          				}
	  				if (is_lexi) {
	    					conf_actual = conf_0;
	    					for (int j = 0; j < n; ++j) {
	      						conf_actual[j][d] = coor_actual[j];
	    					}
	    					//añadir conf_actual y lexi_actual
	    					Gale_conf.push_back(conf_actual);
	    					Gale_lexi.push_back(lexi_actual);
	  				}
        			} while( prev_permutation(coor_actual.begin(), coor_actual.end()));
      			}
    		}
  	}
	//Gale_conf son todas las configuraciones halladas hasta ahora por revisar
	vector <gale_conf> politopes;

	for (unsigned int i=0; i<Gale_conf.size(); ++i){
		//meterlo en una galeconf
		gale_conf G;
		G.vectors=MR(Gale_conf[i]);//Fix this manually to transform ok!!!
		//generar cocircuitos
		set <vector<int> > cocir=generate_cocircuits (G.vectors, n, e);
		bool interior_point=false;
		for (set< vector< int > >::iterator j=cocir.begin(); j!=cocir.end() and !interior_point; ++j){
			int positive=0;
			int negative=0;
			vector<int> v=*j;
			for (int k=0; k<n; ++k){
				if (v[k]>0) ++positive;
				if (v[k]<0) ++negative;
			}
			if ((positive==1 and negative>0) or (positive>0 and negative==1)){//check for interior points
				interior_point=true;
			}
		}
		//from here find if ther are interior points, if there are no need to store the configuration, go to the next one
		if (not interior_point){
			//perform the galetransform
			G.transform=null_space(T(G.vectors));//points of the conf as rows
			//if there is an all 1 column substitute for the first one if not put the first column all ones
			int column=-1;
			for (int j=0; j<n-e and column==-1; ++j){
				bool still_equal=true;
				for(int k=0; k<n-1 and still_equal; ++k){
					if (G.transform[k][j]!= G.transform[k+1][j]){
						still_equal=false;
					}
				}
				if (still_equal){
					column=j;
				}
			}
			if (column==-1){
				for (int j=0; j<n; ++j){
					G.transform[j][0]=1;
				}
			}
			else{
				for (int j=0; j<n; ++j){
					G.transform[j][column]=G.transform[j][0];
					G.transform[j][0]=1;
				}
			}
			//find the cocircuits and compare with the other
			G.cocircuits=generate_cocircuits_trans (G.transform, n, n-e-1);
			bool already_found=false;
			for (unsigned int j=0; j<politopes.size() and !already_found; ++j){
				if (compare_cocircuits(G.cocircuits,politopes[j].cocircuits)){
					already_found=true;
				} 
			}
			if (not already_found) politopes.push_back(G);
			//find the facet lattice and compare with the other
			/*
			else{
				//generate the facet lattice for this gale conf
				G.facets=facet_lattice(G.cocircuits);
				if (coincident_cocircuits.size()==1){
					//generate the facet lattice for the other gale conf
					G.facets=facet_lattice(politopes[coincident_cocircuits[0]].cocircuits); 
					//if there are more than one they have already been generated
				}
				//compare the facet lattices
				for (int j=0; j<coincident_cocircuits.size(); ++j){
					//como esta ordenado puedes comparar directamente
				}
			}
			*/
		}
	}
	return politopes;

}

vector <gale_conf> diference (const  int e, const   int n, const  int m){

	vector <gale_conf> M =galecomplexity (e,n,m);
	vector <gale_conf> previous=galecomplexity(e,n,m-1);
	vector <gale_conf> diference;
	//compare cocircuits and facets (if necessary);
	//provisional solo compara cocircuitos
	for (unsigned int j=0; j<M.size(); ++j){
	  	bool found=false;
	  	for (unsigned int i=0; i<previous.size() and !found; ++i){
	    	if (previous[i].cocircuits==M[j].cocircuits) found=true;
	  	}
	  	if (not found) diference.push_back(M[j]);
	} 
	return diference;
	

}


UserFunctionTemplate4perl("# @category Computations"
			  "# Computes all the possible gale configurations with gale complexity m "
			  "# such that they are not combinatorially equivalent."
			  "# @param e dimension of the gale diagram "
			  "# @param n number of vectors on the gale diagram"
			  "# @param m maximum value of the coordinates"
			  "# @return Array<gale_conf>",
			  "galecomplexity ( int,  int,  int)" );

UserFunctionTemplate4perl("# @category Computations"
			  "# Computes all the possible gale configurations with gale complexity m and m-1 "
			  "# returns only the ones that exist exclusively with complexity m."
			  "# @param e dimension of the gale diagram "
			  "# @param n number of vectors on the gale diagram"
			  "# @param m maximum value of the coordinates"
			  "# @return Array<gale_conf>",
			  "diference ( int ,  int,  int)" );


}}


