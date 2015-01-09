/* Copyright (c) 2014
   Julian Pfeifle & Discrete Geometry class FME/UPC 2014

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

#include "polymake/SparseVector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/ListMatrix.h"
#include "polymake/Array.h"

namespace polymake { namespace polytope {

void print(const Matrix<Integer> &m) {
  for (int i=0; i<m.rows(); ++i) {
    for (int j=0; j<m.cols(); ++j) {
      cout << " " << m[i][j];
    }
    cout << endl;
  }
}

void print(const Matrix<Rational> &m) {
  for (int i=0; i<m.rows(); ++i) {
    for (int j=0; j<m.cols(); ++j) {
      cout << " " << m[i][j];
    }
    cout << endl;
  }
}

void print(const Vector<Rational> &v) {
  for (int i=0; i<v.dim(); ++i) {
      cout << " " << v[i];
  }
  cout << endl;
}




Matrix<Integer> enumerate_configurations(int e, int n, int m) {
  Matrix<Rational> equations(e+1, e*n+1);
  for (int j=0; j<e; ++j) {
    for (int i=0; i<n; ++i) {
      equations(j, 1 + i*e + j) = 1;
    }
  }
  
  // Adding last equation v_00 = m
  equations(e, 1) = -m;
  
  cerr << equations << endl;
  
  ListMatrix<SparseVector<Rational> > inequalities(0, e*n+1);
  
  // Symmetries
  for (int i=0; i<e-1; ++i) { // leave out last row
    SparseVector<Rational> ineq(e*n+1);
    ineq[1 + i] = 1;
    ineq[1 + i + 1] = -1;
    inequalities /= ineq;
  }
  SparseVector<Rational> last(e*n+1);
  last[e] = 1;
  inequalities /= last;
  
  // (Almost) Lexicographic order
  for (int i=0; i<n-1; ++i) { // leave out last row, since it makes no sense
    SparseVector<Rational> ineq(e*n+1);
    ineq[1 + i*e] = 1;
    ineq[1 + (i+1)*e] = -1;
    inequalities /= ineq;
  }
  
  cerr << Matrix<Rational>(inequalities) << endl;
  
  perl::Object P("Polytope");
  P.take("EQUATIONS") << equations;
  P.take("INEQUALITIES") << Matrix<Rational>(inequalities);
  //const Matrix<Integer> C = P.give("LATTICE_POINTS"); // Això no va
  const Matrix<Integer> C = P.CallPolymakeMethod("LATTICE_POINTS"); // S'ha de fer així

  return C;
}


Array<Matrix<Integer> > gale_complexity(int e, int n, int m) {
  Matrix<Integer> mat = enumerate_configurations(e,n,m);
  
  Array<Matrix<Integer> > a(1, mat);

  return a;
}

UserFunctionTemplate4perl("# @category Computations"
                          "# Compute a vector that selects the inclusion-minimal face "
			  "# that contains the vertices indexed by I."
                          "# @param Polytope P the input polytope"
                          "# @param Set I the indices of some vertices of P"
                          "# @return Vector",
			  "gale_complexity(Int, Int, Int)" );

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
