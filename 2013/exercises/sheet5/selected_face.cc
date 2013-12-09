/* Copyright (c) 2013
   ####YOUR_NAME_HERE####

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
#include "polymake/Set.h"
#include "polymake/Rational.h"
#include "polymake/Vector.h"

namespace polymake { namespace polytope {

Set<int> selected_face(perl::Object P, const Vector<Rational>& a)
{
   // your code here
}

UserFunction4perl("# @category Producing from scratch"
                  "# Find the face of P maximizing the linear function a."
                  "# @param Polytope P"
                  "# @param Vector a the linear function to be maximized"
                  "# @return Set",
                  &selected_face, "selected_face(Polytope, Vector)");


} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
