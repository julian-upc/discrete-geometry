/* Copyright (c) 2014
Irene de Parada & Elisabet Burjons Geometry class FME/UPC 2014
This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version: http://www.gnu.org/licenses/gpl.txt.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
--------------------------------------------------------------------------------
*/
#include "polymake/Set.h"
#include "polymake/Vector.h"
#include "polymake/SparseVector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/IncidenceMatrix.h"
#include <set> //compare polymake and c++ sets and decide for one of them

using namespace std;
using namespace polymake;

typedef vector<int> VE; //integer vector
typedef Matrix <Rational> MR;
typedef ListMatrix <SparseVector<Rational>> LMR;

class gale_conf{
	Vector <int> //contains all the vectors one afther the other

	Matrix<Rational> //to be able to calculate kernel

	Matrix<int> //cocircuits

	//vertice set?

	//facet lattice?

	//compare function?
}

void is_new_configuration ( int n, int m, int e, VE gale_conf, /*someclass*/& found_polygons){//looks if a configuration is equivalent to a found polygon and appends it if not


}

MR enumerate_configurations (int e, int n, int m){
	MR equations (e+1,e*n+1);
	for (int j=0; j<e; ++j){
		for (int i=0; i<n; ++i){
			equations(j,i*e+j+1)=1;
		}
	}
	equations(e,0)=-m;
	equations(e,1)=1;

	LMR inequalities (0,e*n+1);// this can be done with matrices and define the size from begining
	for (int i=0; i<e-1; ++i){//leave out last row
		SparseVector<Rational> ineq (e*n+1);
		ineq[i+1]=1;
		ineq[i+2]=-1;
		inequalities/=ineq;
	}
	SparseVector<Rational> ineq (e*n+1);
	ineq[e]=1;
	inequalities/=ineq;

	for (int i=0; i<n-1; ++i){ //yes, really!
		for (int j=0; j<e; ++j){
			SparseVector<Rational> ineq (e*n+1);
			ineq[e*i+j+1]=1;
			ineq[e*(i+1)+j+1]=-1;
			inequalities/=ineq;
		}
	}

	Object P ("Polytope");
	P.take("EQUATIONS")<<equations;
	P.take("INEQUAILITIES")<< MR(inequalities);
	const P.Matrix<Ingeger> C=P.callPolymakeMethod("LATTICE POINTS");
	return C;
}

//falta capçalera!
vector<gale_conf> galecomplexity (int e, int n, int m,  )// I really want to output a set
	//missing input output lines
	const Matrix<Integers> configs= enumerate_configurations(e, n, m);
	for(Entire<Rows<Matrix<Integer>>>::const_iterator rit){

	}
}

/*Vector<Rational> face_selector(perl::Object P, const Set<int>& I)
{
if (I.size() == 0) throw std::runtime_error("The input set of indices must not be empty");
const Matrix<Rational> F = P.give("FACETS");
const IncidenceMatrix<> VIF = P.give("VERTICES_IN_FACETS");
const int ambient_dim (F.cols());
// accumulate the indices of facets that contain the vertices indexed by I
Set<int> facets_containing_I;
for (int i=0; i < VIF.rows(); ++i) {
if (incl(I, VIF[i]) <= 0) {
facets_containing_I += i;
}
}
Vector<Rational> selector(unit_vector<Rational>(ambient_dim, 0));
for (Entire<Set<int> >::const_iterator fit = entire(facets_containing_I); !fit.at_end(); ++fit) {
selector += F[*fit];
}
return selector;
}
UserFunctionTemplate4perl("# @category Computations"
"# Compute a vector that selects the inclusion-minimal face "
"# that contains the vertices indexed by I."
"# @param Polytope P the input polytope"
"# @param Set I the indices of some vertices of P"
"# @return Vector",
"face_selector(Polytope, Set)" );*/

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
