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

perl::Object hypersimplex_slab(int k, int d) {
  if (d < 2)
    throw std::runtime_error("hypersimplex_slab: dimension >= 2 required");
  if (k < 1 || k > d)
    throw std::runtime_error("hypersimplex_slab: 1 <= k <= d required");

  perl::Object p("Polytope<Rational>");
  p.set_description() << "(" << k << "," << d << ")-hypersimplex_slab" << endl;

  p.take("CONE_AMBIENT_DIM") << d + 1;
  p.take("CONE_DIM") << d + 1;

  // we already know the number of vertices
  const int n = Integer::binom(d + 1, k).to_int();
  p.take("N_VERTICES") << n;

  Matrix<Rational> Vertices(n, d + 1);
  Rows< Matrix<Rational> >::iterator v = rows(Vertices).begin();
  Subsets_of_k<sequence> enumerator(range(1, d), k);
  for (Subsets_of_k<sequence>::iterator s = entire(enumerator); !s.at_end(); ++s, ++v) {
    (*v)[0] = 1;
    v->slice(*s).fill(1);
  }
  enumerator = Subsets_of_k<sequence>(range(1, d), k - 1);
  for (Subsets_of_k<sequence>::iterator s = entire(enumerator); !s.at_end(); ++s, ++v) {
    (*v)[0] = 1;
    v->slice(*s).fill(1);
  }
  p.take("VERTICES") << Vertices;
  p.take("LINEALITY_SPACE") << Matrix<Rational>();

  return p;
}

UserFunction4perl("# @category Producing from scratch"
                  "# Produce the hypersimplex slab, that is the the convex hull of all 0/1-vector in R<sup>//d//</sup>"
                  "# with exactly //k// or //k// - 1 1s."
                  "# @param Int k number of 1s"
                  "# @param Int d ambient dimension"
                  "# @return Polytope",
                  &hypersimplex_slab, "hypersimplex_slab");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
