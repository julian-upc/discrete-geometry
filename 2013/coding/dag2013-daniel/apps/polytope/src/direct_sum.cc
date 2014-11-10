/* Copyright (c) 1997-2012
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
   $Project: polymake $$Id: join_polytopes.cc 10553 2012-01-17 18:13:04Z joswig $
*/

#include "polymake/client.h"
#include "polymake/Matrix.h"

namespace polymake { namespace polytope {

template <typename Scalar>
perl::Object join_polytopes(perl::Object p1, perl::Object p2) {
	
   const bool pointed=p1.give("POINTED") && p2.give("POINTED");
   if (!pointed)
      throw std::runtime_error("join_polytopes: input polyhedron not pointed");

   const Matrix<Scalar> v1=p1.give("VERTICES"),
                        v2=p2.give("VERTICES");

   perl::Object p_out(perl::ObjectType::construct<Scalar>("Polytope"));
   const int n1=v1.rows();
   const int n2=v2.rows();
   const Matrix<Scalar> V_out= (v1 | same_element_vector(Scalar(-1),n1) | zero_matrix<Scalar>(n1,v2.cols()-1)) /
                               (ones_vector<Scalar>(n2) | zero_matrix<Scalar>(n2,v1.cols()-1) | v2);

   p_out.set_description()<<"Join of "<<p1.name()<<" and "<<p2.name()<<endl;
   p_out.take("VERTICES") << V_out;
   const Matrix<Scalar> empty;
   p_out.take("LINEALITY_SPACE") << empty;

   return p_out;
}

UserFunctionTemplate4perl("# @category Producing a new polyhedron from others"
                          "# Construct a new polyhedron as the join of two given pointed ones."
                          "# @param Polytope P1"
                          "# @param Polytope P2"
                          "# @return Polytope",
                          "join_polytopes<Scalar>(Polytope<Scalar> Polytope<Scalar>)");
} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
