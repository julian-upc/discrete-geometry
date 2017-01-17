//
//  Gale_Diagram.hpp
//  PolytopeProject
//


#ifndef Gale_Diagram_hpp
#define Gale_Diagram_hpp

#include <stdio.h>
#include "Point.hpp"
#include "Matrix.hpp"
#include "igraph.h"

class Gale_Diagram{
    vector<Point> points;
    Matrix matrix;
    
    vector<int> findMissing(int i, int j, int k,int l) const;
    
    //Does NOT check for degeneracies
    //returns -1 if x lies on the right of the oriented line ab, otherwise +1;
    int orientation(const Point &x,const Point &a,const Point &b) const;
    
    //Checks if onePoint lies in the triangle defined by threePoints
    bool in_triangle(const vector<Point> &onePoint, const vector<Point> &triangle) const;
    
    // returns whether both linesegments (vectors of size 2) intersect.
    bool intersectingLineSegments(const vector<Point> &lineSegment1, const vector<Point> &lineSegment2) const;
    
    bool check_intersecting(const vector<Point> &points) const;
    
    bool isSimplex(int i, int j, int k, int l) const;
public:
    Gale_Diagram(vector<Point> points);
    
    bool is_neighborly() const;
        
    void print() const{
        matrix.print();
    }
    
    //Make the facet-graph.
    igraph_t makeVertexFacetStructure() const;
    
    bool isSimplicial() const;
    
    void test() const;
    
    bool isIsomorphic(const Gale_Diagram &gale) const;
    
    void writeGraph(string filename);
    
    Matrix galeToPolytope() const;
};

#endif /* Gale_Diagram_hpp */
