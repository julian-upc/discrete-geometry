/* Copyright (c) 2013
   Anna Somoza

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
perl::Object direct_sum(perl::Object p, perl::Object q)
{
   const bool pointed=p.give("POINTED") && q.give("POINTED");
   if (!pointed)
      throw std::runtime_error("direct_sum: input polyhedron not pointed");

   const Matrix<Scalar> v = p.give("VERTICES"), w=q.give("VERTICES");

   perl::Object p_out(perl::ObjectType::construct<Scalar>("Polytope"));
   const int n = v.rows(), m = w.rows();
   const Matrix<Scalar> v_out = (v | zero_matrix<Scalar>(n, w.cols()-1)) /
      (ones_vector<Scalar>(m) | zero_matrix<Scalar>(m,v.cols()-1) | w.minor(All, range(1, w.cols()-1)));
   p_out.set_description() << "Direct Sum of " << p.name() << " and " << q.name() << endl;
   p_out.take("VERTICES") << v_out;
   const Matrix<Scalar> empty;
   p_out.take("LINEALITY_SPACE") << empty;

   return p_out;
}

UserFunctionTemplate4perl("# @category Producing a new polyhedron from others"
                          "# Construct a new polyhedron as the direct sum of two pointed ones."
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
