/* Copyright (c) 2013
   Alex Alvarez Ruiz, Universitat Politecnica de Catalunya

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

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/PowerSet.h"

namespace polymake { namespace polytope {

perl::Object sigma(int k, int d) {
  if (d < 2)
    throw std::runtime_error("sigma: dimension >= 2 required");
  if (k < 1 || k >= d)
    throw std::runtime_error("sigma: 1 <= k < d required");

  perl::Object p("Polytope<Rational>");
  p.set_description() << "(" << k << "," << d << ")-hypersimplex_slab" << endl;

  p.take("CONE_AMBIENT_DIM") << d + 1;
  p.take("CONE_DIM") << d;

  // we already know the number of vertices
  const int n = 1 + k*(d - k);
  p.take("N_VERTICES") << n;

  Matrix<Rational> Vertices(n, d + 1);
  Rows< Matrix<Rational> >::iterator v = rows(Vertices).begin();
  
  //The vertex of the simplex that is intersected
  for (int i = 0; i <= k; ++i) (*v)[i] = 1;
  ++v;
  
  //Compute the intersections with the edges
  for (int i = 0; i < k; ++i) {
    for (int j = k + 1; j <= d; ++j) {
      (*v)[0] = 1;
      const Rational lambda(k - i, j - i);
      for (int t = 0; t < d; ++t) {
        const int u = (i > t ? 1 : 0), w = (j > t ? 1 : 0);
        (*v)[t + 1] = u + lambda*(w - u); 
      }
      ++v;
    }
  }
  p.take("VERTICES") << Vertices;
  p.take("LINEALITY_SPACE") << Matrix<Rational>();

  return p;
}

UserFunction4perl("# @category Producing from scratch"
                  "# Produce the Sigma slice, that is the convex hull of the intersection of the hyperplane \\sum x_i = k"
                  "# with the polytope defined by conv(0, e_1, e_1 + e_2, ..., e_1 + ... + e_d) in R<sup>//d//</sup>"
                  "# @param Int k number of 1s"
                  "# @param Int d ambient dimension"
                  "# @return Polytope",
                  &sigma, "sigma");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
