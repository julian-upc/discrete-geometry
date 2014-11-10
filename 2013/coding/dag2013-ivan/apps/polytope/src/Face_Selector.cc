/* Copyright (c) 2013
Julian Pfeifle

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

#include "polymake/client.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Set.h"

namespace polymake { namespace polytope {

Vector<Rational> face_selector(perl::Object P, const Set<int>& G)
{
    const IncidenceMatrix<> VIF = P.give("VERTICES_IN_FACETS");
    const Matrix<Rational> facets = P.give("FACETS");
    const int ambient_dim = facets.cols();

    // sum the normal vectors corresponding to the facets containing G
    Vector<Rational> selector(unit_vector<Rational>(ambient_dim, 0));
    int i(0);
    for (Entire<Rows<IncidenceMatrix<> > >::const_iterator rit = entire(rows(VIF)); !rit.at_end(); ++rit, ++i) {
        const Set<int> F(*rit);
        if (incl(G,F) == -1)
         selector += facets[i];
    }

    return selector;
}

UserFunction4perl("# @category Producing from scratch"
                  "# Produce a linear objective function that selects the smallest face containing G."
                  "# @param Polytope P "
                  "# @param Set G index set"
                  "# @return Vector",
                  &face_selector, "face_selector(Polytope, Set)");


} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End
