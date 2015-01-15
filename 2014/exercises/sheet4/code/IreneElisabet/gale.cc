//Darrera edició: 14 de gener de 2015

#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <math.h> 
#include <algorithm> 
#include "polymake/Rational.h"
#include "polymake/Vector.h"
#include "polymake/Set.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/linalg.h"


using namespace std;

namespace polymake{ namespace polytope{

typedef Matrix <Rational> MR;
typedef vector< vector< vector<int> > > Tri_vector;
typedef vector< vector<int> >  Matriz;



struct gale_conf{
	Matrix<Rational> vectors; //to be able to calculate kernel breaks the configuration into vectors

	Matrix <Rational> transform; //stores the transformation of the gale configuration

	set<vector<int> > cocircuits; //cocircuits of the transform

	set< set<int> > facets; //stores the facet lattice cocircuits of the transform
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

void write_matrix (const MR & A) {
  for (int i = 0; i < A.rows (); ++i) {
    for (int j = 0; j < A.cols(); ++j) {
      cout << A[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
}

void write_vector ( const vector <int>& z) {
  for (unsigned int i = 0; i < z.size(); ++i) {
    cout << z[i] << " ";
  }
  cout << endl;
}

void write_stdintmatrix (const Matriz & M){
    for (unsigned int i = 0; i < M.size(); ++i) {
    write_vector(M[i]);
  }
  cout<<endl;
}

vector <Rational> hiperplane ( const MR & M){
	MR x=null_space (M);
	vector<Rational> x2 (x.cols());
	for (int i=0; i<x.cols(); ++i){
		x2[i]=x[0][i];
	}
	return x2;
}



set < vector<int > > generate_cocircuits ( const Matrix <Rational> & M, const  int n, const  int e){
	Matriz c = combinations (n, e-1);//han de passar pel zero
	set<vector<int> > cocircuits;
	for (unsigned int i=0; i<c.size(); ++i){
		MR A (e-1,e);
		for (int j=0; j<e-1; ++j){
			for(int k=0; k<e; ++k){
				A[j][k]=M[c[i][j]][k];
			}
		}
		if (rank(A)==e-1){//make sure we can build a hiperplane
			vector< Rational> x= hiperplane(A);//hiperplane
			vector<int> y(n);
			for (int j=0; j<n; ++j){
				Rational a=0;
				for (int k=0; k<e; ++k){
					a=a+x[k]*M[j][k];
				}
				if (a>0) y[j]=1;
				else if (a==0) y[j]=0;
				else y[j]=-1;
			}
			cocircuits.insert(y);
		}
	}
	return cocircuits;
}


set<vector< int > > generate_cocircuits_trans ( const Matrix <Rational> & M, const  int n, const  int d){
	Matriz c = combinations (n, d);
	set<vector < int > > cocircuits;
	for (unsigned int i=0; i<c.size(); ++i){
		MR A (d,d+1);
		for (int j=0; j<d; ++j){
			for(int k=0; k<d+1; ++k){
				A[j][k]=M[c[i][j]][k];
			}
		}
		if (rank(A)==d){
			vector< Rational> x= hiperplane(A);
			vector<int> y(n);
			for (int j=0; j<n; ++j){
				Rational a=0;
				for (int k=0; k<d+1; ++k){
					a=a+x[k]*M[j][k];
				}
				if (a>0) y[j]=1;
				else if (a==0) y[j]=0;
				else y[j]=-1;
			}
			cocircuits.insert(y);
		}
	}
	return cocircuits;
}



set <set <int> > facet_lattice (const set<vector< int> > & cocircuits){
	set <set <int> > facets;
	for (set<vector< int > >::iterator it= cocircuits.begin(); it !=cocircuits.end(); ++it){
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
			set<int> facet;
			for (unsigned int i=0; i<v.size(); ++i){//with this all the facets are positive
				if (v[i]==0) facet.insert(i);
			}
			facets.insert(facet);
		}
	}
	return facets;
}



Matriz permuta (int n){
	vector <int> v (n);
	Matriz p;
	for (int i=0; i<n; ++i){
		v[i]=i;
	}
	do {
		p.push_back(v);
	} while ( next_permutation(v.begin(),v.end()) );
	return p;
}

set <set <int> > perm (const  set< set< int> > & f, const vector< int> & p){
	set <set< int> > f2;
	for (set < set< int> >::iterator it= f.begin(); it!= f.end(); ++it){
		set <int> facet=*it;
		set <int> new_facet;
		for (set <int>:: iterator fit=facet.begin(); fit!= facet.end(); ++fit){
			int a=*fit;
			new_facet.insert(p[a]);
		}
		f2.insert(new_facet);
	}
	return f2;
}



bool compare_facets (const set< set <int> > & f1, const set< set <int> > & f2, const int & n){
	if (f1==f2) return true;
	Matriz p=permuta (n);
	for (unsigned int i=1; i<p.size(); ++i){//the first element of p is the non permutated one
		set< set< int> > changed= perm (f2, p[i]);
		if (f1==changed) return true;
	}
	return false;
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
		Coor1D.insert(coor_actual);
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
	int full_dim=0;
	int intpoint=0;
	for (unsigned int i=0; i<Gale_conf.size(); ++i){
		//meterlo en una galeconf
		gale_conf G;
		G.vectors=MR(n,e);//Fix this manually to transform ok!!!
		for (int k=0; k<e; ++k){
			for (int j=0; j<n; ++j){
				G.vectors[j][k]=Rational(Gale_conf[i][j][k]);
			}
		}
		if (rank(G.vectors)==e){
			++full_dim;
			//generar cocircuitos
			set <vector<int > > cocir=generate_cocircuits (G.vectors, n, e);
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
				++intpoint;
				//perform the galetransform
				//add 1's row to T(G.vectors);
				MR M (e+1, n);
				for (int j=0; j<n; ++j){//1st row
					M[0][j]=1;
				}
				for (int j=0; j<n; ++j){//copy the transposed matrix
					for (int k=0; k<e; ++k){
						M[k+1][j]=G.vectors[j][k];
					}
				}
				//write_matrix(M);
				MR A=T(null_space(M));
				//cout << "transformed,  not interior point" << endl;
				//add a 1's column
				int d=n-e-1;
				MR A2 (n, d+1);
				for (int j=0; j<n; ++j){
					A2[j][0]=1;
				}
				for (int j=0; j<n; ++j){
					for (int k=0; k<d; ++k){
						A2[j][k+1]=A[j][k];
					}
				}
				G.transform=A2;
				
				//find the cocircuits and compare with the other
				G.cocircuits=generate_cocircuits_trans (G.transform, n, d);
				//find the facet lattice and compare with the other
				//alternative for facets
			
				G.facets=facet_lattice(G.cocircuits);
				bool already_found_facets=false;
				for (unsigned int j=0; j<politopes.size() and !already_found_facets; ++j){
					if (compare_facets(G.facets,politopes[j].facets, n)){
						already_found_facets=true;
					} 
				}
				if (not already_found_facets){
					politopes.push_back(G);
					/*
					write_matrix (G.transform);
					for (set<set<int> >::iterator it=G.facets.begin(); it!=G.facets.end(); ++it){
						set<int> v=*it;
						for (set<int>::iterator fit=v.begin(); fit!=v.end(); ++fit){
							int punt=*fit;
							cout<<punt<<" ";
						}
						cout<<endl;
					}
					cout<<endl;
					*/
				}
			}
			
		}
		
	}
	/*
	cout<<"tengo "<<Gale_conf.size()<<"posibles configuraciones"<<endl;
	cout<<" de ellas "<<full_dim<<" son full dimensional"<<endl;
	cout<<" de ellas "<<intpoint<<" no tienen puntos interiores"<<endl;
	cout<<" en total "<<politopes.size()<<" configuraciones son validas"<<endl;
	*/
	return politopes;

}

Vector<Matrix<Rational > > convert (const  vector<gale_conf> & M ){
	Vector<Matrix<Rational> > v;
	v.clear();
	for (unsigned int i=0; i<M.size(); ++i){
		v|=M[i].vectors;//appends A to the vector
	}
	return v;
}

Vector<Matrix<Rational> > gale_complexity (const  int e, const   int n, const  int m){
	vector<gale_conf> M =galecomplexity (e,n,m);
	Vector< Matrix <Rational> > politopes=convert (M);
	return politopes;
}


Vector<Matrix< Rational > > diference (const  int e, const   int n, const  int m){

	vector <gale_conf> M =galecomplexity (e,n,m);
	vector <gale_conf> previous=galecomplexity(e,n,m-1);
	vector <gale_conf> diference;
	//compare cocircuits and facets (if necessary);
	//provisional solo compara cocircuitos
	for (unsigned int j=0; j<M.size(); ++j){

		//alternative for facets
		
		bool found_facets=false;
	  	for (unsigned int i=0; i<previous.size() and !found_facets; ++i){
		      if (compare_facets( previous[i].facets, M[j].facets, n)) found_facets=true;
	  	}
	  	if (not found_facets) diference.push_back(M[j]);
	  	
	} 
	
	return convert(diference);

}


UserFunctionTemplate4perl("# @category Computations"
			  "# Computes all the possible gale configurations with gale complexity m "
			  "# such that they are not combinatorially equivalent."
			  "# @param e dimension of the gale diagram "
			  "# @param n number of vectors on the gale diagram"
			  "# @param m maximum value of the coordinates"
			  "# @return Vector<Matrix<Rational> >",
			  "gale_complexity ($ $ $)" );

UserFunctionTemplate4perl("# @category Computations"
			  "# Computes all the possible gale configurations with gale complexity m and m-1 "
			  "# returns only the ones that exist exclusively with complexity m."
			  "# @param e dimension of the gale diagram "
			  "# @param n number of vectors on the gale diagram"
			  "# @param m maximum value of the coordinates"
			  "# @return Vector<Matrix<Rational> >",
			  "diference ($ $ $)" );


}}

