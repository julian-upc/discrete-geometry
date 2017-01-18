#include <iostream>
#include <vector>
#include <cmath>

using namespace std;
typedef double num;
//This represents  a point in the plane
class Point {
  num x;
  num y;
public:
  Point();
  Point(double _x, double _y);
  num operator[](int index) const {
    if(index == 0) return x;
    else if(index == 1) return y;
    else return 133;
  }
  num getX() const {return x;}
  num getY() const {return y;}
} ;
Point::Point(){
  x = 1000;
  y = 1000;
}
Point::Point(num _x, num _y){
  this-> x = _x;
  this-> y = _y;
}

//Represents a single "point configuration". For each point, it can be negative or positive, in position 0 or in position 1
class PointConfig {
  vector<Point> points;
public:
  bool positive;
  bool consider;
  int index_point;
  PointConfig();
  PointConfig(Point p1, Point p2);
  const Point& currentPoint() const;
  const Point& operator[](int index) const {
    return this->points[index];
  }
  void print();

} ;
PointConfig::PointConfig(){
  PointConfig(Point(),Point());
}
PointConfig::PointConfig (Point p1, Point p2){
  this->points = vector<Point>(2);
  this->points[0] = p1;
  this->points[1] = p2;
  this->positive = true;
  this->consider = true;
  this->index_point = 0;
}
 const Point& PointConfig::currentPoint() const{
  return this->points[this->index_point];
}

void PointConfig::print(){
  Point p = this->currentPoint();
  cout << " ({" << p[0] << "," << p[1] << "}, " << ((this->positive)?"+":"-") << ") ";
}
//Helper class to compute intersections when calculating whehter the convex hull of positive points
//intersects the convex hull of negative points
class LineSegment {
  Point a;
  Point b;
public:
  LineSegment(Point _a, Point _b);
  const Point& operator[](int index) const {
    if(index == 0) return this->a;
    else return this->b;
  }
};

LineSegment::LineSegment(Point _a, Point _b) : a(Point(_a[0],_a[1])), b(Point(_b[0],_b[1])){}

//Computes the signed area of p,q,r and returns the orientation from that as seen in class
num orientation2D(const Point& p, const Point& q, const Point& r)
{
  num qpx = q[0]-p[0];
  num rpy = r[1]-p[1];
  num rpx = r[0]-p[0];
  num qpy = q[1]-p[1];

  num result = qpx*rpy - rpx*qpy;
  if(result > 0) return 1;
  else if(result < 0) return -1;
  else return 0;
}

enum Intersection {EMPTY, INTERIOR, SAME, ENDPOINT_OVERLAP, OVERLAP,
                   PARTIAL_OVERLAP, ENDPOINT_ALIGN, ENDPOINT, ENDPOINT_INTERIOR};
//This is a function adapted from Vera's FIB course, originally intended to classify types
//of intersections between line segments, but in our particular case we are only interested
//in whether they intersect or not
Intersection classifyIntersection(const LineSegment& s, const LineSegment& t)
{
    num t0wrtS = orientation2D(s[0],s[1],t[0]);
    num t1wrtS = orientation2D(s[0],s[1],t[1]);
    num s0wrtT = orientation2D(t[0],t[1],s[0]);
    num s1wrtT = orientation2D(t[0],t[1],s[1]);
    Intersection result = EMPTY;
    //If the two segments are aligned
    if(t0wrtS == 0 and t1wrtS == 0 and s0wrtT == 0 and s1wrtT == 0){
      cout << "Shouldn't be here" << endl;
      result = SAME;
      /*LineSegment leftmost = s;
        LineSegment rightmost = t;
        if(leftmost[1] < leftmost[0]) leftmost = leftmost.opposite();
        if(rightmost[1] < rightmost[0]) rightmost = rightmost.opposite();
        if(t.min() < s.min()){
            LineSegment aux = leftmost;
            leftmost = rightmost;
            rightmost = aux;
        }
        if (leftmost[0] == rightmost[0] and leftmost[1] == rightmost[1]){
            result = SAME;
        }else if(leftmost[0] == rightmost[0] or leftmost[1] == rightmost[1]){
            result = ENDPOINT_OVERLAP;
        }else if(leftmost[1] > rightmost[0]){
            if(rightmost[1] <= leftmost[1]) {
               result = OVERLAP; 
            }else {
               result = PARTIAL_OVERLAP;
            } 
        }else if(leftmost[1] == rightmost[0]) {
             result = ENDPOINT_ALIGN;
             }*/
    //Although now we know they are not aligned, they could still share an endpoint
    }else if((t0wrtS == 0 or t1wrtS == 0) and (s0wrtT == 0 or s1wrtT == 0) 
        and (t0wrtS != 0 or t1wrtS != 0 or  s0wrtT != 0 or s1wrtT != 0)){
        result = ENDPOINT;
    //Case where both lines intersect (then we must see if the segments intersect)
    }else if(t0wrtS != 0 and t1wrtS != 0 and  s0wrtT != 0 and s1wrtT != 0){
        if(t0wrtS != t1wrtS and s0wrtT != s1wrtT) {
            result = INTERIOR;
        }
    //Case where one endpoint is intersecting with interior points of the other
    //segment
    }else if((t0wrtS == 0 or t1wrtS == 0) and s0wrtT != s1wrtT){
            result = ENDPOINT_INTERIOR;
    }else if((s0wrtT == 0 or s1wrtT == 0) and t0wrtS != t1wrtS){
            result = ENDPOINT_INTERIOR;
    }
    return result;
}


