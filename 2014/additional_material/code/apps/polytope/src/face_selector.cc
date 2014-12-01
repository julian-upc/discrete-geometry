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

#include "polymake/Set.h"
#include "polymake/Vector.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"

namespace polymake { namespace polytope {

Vector<Rational> face_selector(perl::Object P, const Set<int>& I)
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
			  "face_selector(Polytope, Set)" );

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
