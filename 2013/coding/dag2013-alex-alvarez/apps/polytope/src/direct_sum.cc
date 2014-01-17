/* Copyright (c) 2013
   Alex Alvarez Ruiz, Universitat Politecnica de Catalunya.

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
#include "polymake/Matrix.h"

namespace polymake { namespace polytope {

template <typename Scalar>
perl::Object direct_sum(perl::Object P, perl::Object Q) {
   const bool pointed = P.give("POINTED") and Q.give("POINTED");
   if (not pointed)
      throw std::runtime_error("direct_sum: input polyhedron not pointed");

   const Matrix<Scalar> v1 = P.give("VERTICES"), v2 = Q.give("VERTICES");

   perl::Object p_out(perl::ObjectType::construct<Scalar>("Polytope"));
   const int n1 = v1.rows(), n2 = v2.rows();
   const Matrix<Scalar> V_out = (v1 | zero_matrix<Scalar>(n1, v2.cols() - 1)) /
      (ones_vector<Scalar>(n2) | zero_matrix<Scalar>(n2, v1.cols() - 1) | v2.minor(All, range(1, v2.cols() - 1)));
   
   p_out.set_description() << "Direct Sum of " << P.name() << " and " << Q.name() << endl;
   p_out.take("VERTICES") << V_out;
   const Matrix<Scalar> empty;
   p_out.take("LINEALITY_SPACE") << empty;

   return p_out;
}

UserFunctionTemplate4perl("# @category Producing a new polyhedron from others"
                          "# Construct a new polyhedron as the direct sum of two given pointed ones."
                          "# @param Polytope P"
                          "# @param Polytope Q"
                          "# @return Polytope",
                          "direct_sum<Scalar>(Polytope<Scalar> Polytope<Scalar>)");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
