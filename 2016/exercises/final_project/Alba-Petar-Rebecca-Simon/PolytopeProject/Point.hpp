//
//  Point.hpp
//  PolytopeProject
//


#ifndef Point_hpp
#define Point_hpp

#include <stdio.h>
#include "iostream"
#include <vector>
#include <gmpxx.h>

using namespace std;


class Point{
    bool sign;
    vector<mpq_class> coordinates;
public:
    Point(vector<mpq_class> coord,bool sign);
    
    Point(vector<double> coord, bool sign);
    
    bool get_sign() const;
    void set_sign(bool s);
    vector<mpq_class> get_coordinates() const;
    
    
    //Note that we will present a vector with the Point class.
    Point point_to_vector(Point point);
    
    //returns true if all coordinates have the same sign.
    bool one_sign() const;
    
    //Substract the coordinates of point b from this point and return a new point
    Point minus (const Point &b) const;
    
    //Multiply every coordinate of this point by a scalar and return this in a different point.
    Point multiply(const mpq_class &scalar) const;
    
    void print() const;
    
    void printLatex() const;
    
};


#endif /* Point_hpp */
