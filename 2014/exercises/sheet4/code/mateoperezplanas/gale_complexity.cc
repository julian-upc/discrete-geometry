
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

//#include "polymake/SparseVector.h"
//#include "polymake/ListMatrix.h"
//#include "polymake/SparseMatrix.h"
#include "polymake/Rational.h"
#include "polymake/Matrix.h"
#include "polymake/Array.h"
#include "polymake/Map.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/permutations.h"

namespace polymake { namespace polytope {

class Configuration {
  public:
    int e;
    int n;
    int d;
    Vector<Rational> v;
    Matrix<Rational> m;
    Matrix<Rational> m_dual;
    Matrix<Rational> cocircuits; // Cocircuits of G are circuits of P
    Matrix<Rational> circuits; // Circuits of G are cocircuits of P
    
    Configuration(int _e, int _n, const Vector<Rational>& _v) {
      e = _e;
      n = _n;
      d = n - e -1;
      
      v = _v;
      
      m = Matrix<Rational>(e+1, n);
      for (int j=0; j<n; ++j) {
        m[0][j] = 1;
        for (int i=1; i<e+1; ++i) {
          m[i][j] = v[j*e + i];
        }
      }
      
      m_dual = ones_matrix<Rational>(1, n); // We add ones to use homogeneous coordinates
      m_dual /= null_space(m);
      
      m = m.minor(sequence(1, e), sequence(0, n));
      
      gen_cocircuits();
    }
    
    void gen_cocircuits() {
      Set<int> s = sequence(0, n);
      int k = e-1;
      Subsets_of_k<const Set<int>&> subsets_of_k = all_subsets_of_k(s, k);

      int num_cocircuits = subsets_of_k.size();
      cocircuits = Matrix<Rational>(num_cocircuits, n);

      Matrix<Rational> Tm = T(m);
      
      int l = 0;
      
      // Iterate over all subsets of {0, ..., n-1} of size k
      for (Entire<Subsets_of_k<const Set<int>&> >::const_iterator subset = entire(subsets_of_k); !subset.at_end(); ++subset) {
        Set<int> comb = *subset;
        
        // Generate a matrix to perform determinants in order to know the orientation
        Matrix<Rational> orientation (k+1, k+1);
        int i=0;
        for (Entire<Set<int> >::const_iterator index = entire(comb); !index.at_end(); ++index) {
          orientation[i] = Tm[*index];
          ++i;
        }
        
        // Save the cocircuit. Points in comb have sign 0. For the rest we have to take the sign of the determinant.
        for (int j=0; j<n; ++j) {
          if (comb.contains(j)) {
            cocircuits[l][j] = 0;
          } else {
            orientation[k] = Tm[j];
            Rational d = det(orientation);
            if (d < 0) {
              cocircuits[l][j] = -1;
            } else if (d > 0){
              cocircuits[l][j] = 1;
            } else {
              cocircuits[l][j] = 0;
            }
          }
        }
        
        ++l;
      }
    }
    
    void gen_circuits() {
      Set<int> s = sequence(0, n);
      int k = d;
      Subsets_of_k<const Set<int>&> subsets_of_k = all_subsets_of_k(s, k);

      int num_circuits = subsets_of_k.size();
      circuits = Matrix<Rational>(num_circuits, n);
      
      // We generate circuits as cocircuits of the dual
      Matrix<Rational> Tm_dual = T(m_dual);
      
      int l = 0;
      
      // Iterate over all subsets of {0, ..., n-1} of size k
      for (Entire<Subsets_of_k<const Set<int>&> >::const_iterator subset = entire(subsets_of_k); !subset.at_end(); ++subset) {
        Set<int> comb = *subset;
        
        // Generate a matrix to perform determinants in order to know the orientation
        Matrix<Rational> orientation (k+1, k+1);
        int i=0;
        for (Entire<Set<int> >::const_iterator index = entire(comb); !index.at_end(); ++index) {
          orientation[i] = Tm_dual[*index];
          ++i;
        }
        
        // Save the circuit. Points in comb have sign 0. For the rest we have to take the sign of the determinant.
        for (int j=0; j<n; ++j) {
          if (comb.contains(j)) {
            circuits[l][j] = 0;
          } else {
            orientation[k] = Tm_dual[j];
            Rational d = det(orientation);
            if (d < 0) {
              circuits[l][j] = -1;
            } else if (d > 0){
              circuits[l][j] = 1;
            } else {
              circuits[l][j] = 0;
            }
          }
        }
        
        ++l;
      }
    }
};

class Polytope {
  public:
    int d;
    int n;
    Set<Set<int> > facets;
    int n_facets;
    int max_vertices_on_facet;
    
    Polytope() {}
   
    Polytope(int _d, int _n, const Matrix<Rational>& cocircuits) {
      d = _d;
      n = _n;
      max_vertices_on_facet = 0;
      for (int i=0; i<cocircuits.rows(); ++i) {
        Rational val = 0;
        bool is_facet = true;
        Set<int> facet;
        
        for (int j=0; j<cocircuits.cols(); ++j) {
          if (cocircuits[i][j] == 0) {
            facet.insert(j);
          }
        
          if (val == 0) {
            val = cocircuits[i][j];
          } else if (cocircuits[i][j] != 0 && val != cocircuits[i][j]) {
            is_facet = false;
          }
        }
        
        if (is_facet) {
          int vertices_on_facet = facet.size();
          if (vertices_on_facet > max_vertices_on_facet) {
            max_vertices_on_facet = vertices_on_facet;
          }
          facets.insert(facet);
        }
      }
      
      n_facets = facets.size();
    }
    
    Polytope(const Polytope& p, const std::vector<int> v) {
      d = p.d;
      n = p.n;
      n_facets = p.n_facets;
      max_vertices_on_facet = p.max_vertices_on_facet;
      
      for (Entire<Set<Set<int> > >::const_iterator fit = entire(p.facets); !fit.at_end(); ++fit) {
        Set<int> p_facet = *fit;
        Set<int> facet;
        for (Entire<Set<int> >::const_iterator vit = entire(p_facet); !vit.at_end(); ++vit) {
          facet.insert(v[*vit]);
        }
        facets.insert(facet);
      }
    }
    
    bool operator<(Polytope other) const {
        if (n_facets < other.n_facets) {
          return true;
        }
        
        if (n_facets > other.n_facets) {
          return false;
        }
        
        if (max_vertices_on_facet < other.max_vertices_on_facet) {
          return true;
        }
        
        if (max_vertices_on_facet > other.max_vertices_on_facet) {
          return false;
        }
        
        return false;
    }
    
    bool operator==(Polytope other) const {
      //if (n_facets == other.n_facets && max_vertices_on_facet == other.)
      return false;
    }

};

Matrix<Rational> enumerate_configurations(int e, int n, int m) {
  Matrix<Rational> equations(e+1, e*n+1);
  for (int j=0; j<e; ++j) {
    for (int i=0; i<n; ++i) {
      equations(j, 1 + i*e + j) = 1;
    }
  }
  
  // Adding last equation v_11 = m
  equations(e, 0) = -m;
  equations(e, 1) = 1;
  
  Matrix<Rational> inequalities(2*e*n + n + e - 1, e*n+1);
  
  // Symmetries
  int index = 0;
  for (int i=0; i<e-1; ++i) { // leave out last row
    inequalities[index + i][1 + i] = 1;
    inequalities[index + i][1 + i + 1] = -1;
  }
  inequalities[e-1][e] = 1;
  index += e;
  
  // (Almost) Lexicographic order
  for (int i=0; i<n-1; ++i) { // leave out last row, since it makes no sense
    inequalities[index + i][1 + i*e] = 1;
    inequalities[index + i][1 + (i+1)*e] = -1;
  }
  index += n-1;
  
  // Bounding box inequalities
  for (int i=1; i<e*n+1; ++i) {
    inequalities[index + 2*(i - 1)][0] = m;
    inequalities[index + 2*(i - 1)][i] = 1;
    inequalities[index + 2*i - 1][0] = m;
    inequalities[index + 2*i - 1][i] = -1;
  }
  
  perl::Object P("Polytope");
  P.take("EQUATIONS") << equations;
  P.take("INEQUALITIES") << inequalities;
  
  const Matrix<Rational> C = P.CallPolymakeMethod("LATTICE_POINTS"); // S'ha de fer així
  
  return C;
}

bool check_lexicographically_sorted(int e, const Vector<Rational>& v) {
  int n = (v.size()-1)/e;

  int prev = 1;
  for (int i=1; i<n; i++) {
    int next = 1 + e*i;
    
    for (int j=0; j<e; ++j) {
      if (v[prev + j] > v[next + j]) {
        break;
      }
      
      if (v[prev + j] < v[next + j]) {
        return false;
      }
    }
    
    prev = next;
  }
  
  return true;
}

bool check_internal_points(const Matrix<Rational>& m) {
  for (int i=0; i<m.rows(); ++i) {
    int pos = 0;
    int neg = 0;
    for (int j=0; j<m.cols(); ++j) {
      if (m[i][j] > 0) pos++;
      if (m[i][j] < 0) neg++;
    }
    
    if (pos < 2 || neg < 2) {
      return true;
    }
  }
  
  return false;
}

bool check_equivalent(const Polytope& p, const Polytope& q) {
  if (p.facets == q.facets) {
    return true;
  }
  
  if (p.d != q.d || p.n != q.n) {
    return false;
  }
  
  AllPermutations<> all_perm = all_permutations(p.n);
  for (Entire<AllPermutations<> >::const_iterator it = entire(all_perm); !it.at_end(); ++it) {
    std::vector<int> v = *it;
    
    Polytope r = Polytope(q, v);
    
    if (p.facets == r.facets) {
      return true;
    }
  }
  
  return false;
}

Array<Matrix<Rational> > gale_complexity(int e, int n, int m) {
  Matrix<Rational> configs = enumerate_configurations(e,n,m);
  
  Map<int, Map<int, Set<Polytope> > > polytopes;
  
  for (int k=0; k<configs.rows(); ++k) {
    if (check_lexicographically_sorted(e, configs[k])) {
      Configuration conf(e, n, configs[k]);
      
      // Internal points filter
      if (!check_internal_points(conf.cocircuits)) {
        conf.gen_circuits();
         
        Polytope p(conf.d, n, conf.circuits);
        
        bool is_new = true;
        for (Entire<Set<Polytope> >::const_iterator it = entire(polytopes[p.n_facets][p.max_vertices_on_facet]); !it.at_end() && is_new; ++it) {
          Polytope current = *it;
          if (check_equivalent(p, current)) {
            is_new = false;
          }
        }
        
        if (is_new) {
          polytopes[p.n_facets][p.max_vertices_on_facet].insert(p);
        }
      }
    }
  }
  
  int total = 0;
  for (Entire<Map<int, Map<int, Set<Polytope> > > >::const_iterator mitf = entire(polytopes); !mitf.at_end(); ++mitf) {
    for (Entire<Map<int, Set<Polytope> > >::const_iterator mitv = entire(mitf->second); !mitv.at_end(); ++mitv) {
      cout << mitf->first << " facets, " <<  mitv->first << " max vertices on facet:" <<  endl;
      for (Entire<Set<Polytope> >::const_iterator sit = entire(mitv->second); !sit.at_end(); ++sit) {
        cout << (*sit).facets << endl;
        total++;
      }
    }
  }
  cout << "Total: " << total << endl;
  
  return Array<Matrix<Rational> > (0);
}

UserFunction4perl("# @category Computations"
                  "# Computes the lattice points of all polytopes with Gale complexity m."
                  "# @param Int e Dimension of the Gale diagram vector space"
                  "# @param Int n Number of vectors of the Gale diagram (i.e. number of vertices of the polytopes)"
                  "# @param Int m Gale complexity of the polytopes"
                  "# @return Array<Matrix<Rational> > Array of all the lattice points for all the polytopes",
                  &gale_complexity,
                  "gale_complexity($$$)" );

} }

// Local Variables:
// mode:C++
// c-basic-offset:3
// indent-tabs-mode:nil
// End:
