#include <iostream>
#include <vector>
#include <cmath>
//#include "auxiliar.cpp"
#include "diedral.cpp"

using namespace std;
//Test whether the convex hulls of the positives and that of the negatives intersect
bool configIntersects(vector<PointConfig>& c){
  vector<PointConfig> positives = vector<PointConfig>();
  vector<PointConfig> negatives = vector<PointConfig>();

  for(unsigned long i = 0; i < c.size(); ++i){
    if(c[i].consider){
      if(c[i].positive) positives.push_back(c[i]);
      else negatives.push_back(c[i]);
    }
  }
  vector<LineSegment> posSeg = vector<LineSegment>();
  vector<LineSegment> negSeg = vector<LineSegment>();
  //Build line segments for the positive points
  for(unsigned long i = 0; i < positives.size(); ++i) {
    for(unsigned long j = i+1; j < positives.size(); ++j){
      posSeg.push_back(LineSegment(positives[i].currentPoint(), positives[j].currentPoint()));
    }
  }
  //Build line segments for the negative points
  for(unsigned long i = 0; i < negatives.size(); ++i) {
    for(unsigned long j = i+1; j < negatives.size(); ++j){
      negSeg.push_back(LineSegment(negatives[i].currentPoint(), negatives[j].currentPoint()));
    }
  }
  //Check if any of the previous line segments, intersect. This is not the smartest way,
  //it takes quadratic time, but the size of the convex hulls is fixed for our problem,
  //since the input is the number of potential configurations, so we could regard this operation as
  //constant time
  bool result = false;
  for(unsigned long i = 0; i < posSeg.size(); ++i) {
    for(unsigned long j = 0; j < negSeg.size(); ++j){
      Intersection aux = classifyIntersection(posSeg[i], negSeg[j]);
      if(aux != EMPTY) return true;
    }
  }
  return result;
}
//Helper function to print the configuration
  void printConfig(vector<PointConfig>& c){
    for(unsigned long i = 0; i < c.size(); ++i){
      c[i].print();
    }
    cout << endl;
  }
//Self explanatory, checks if the current configuration could be that of a Neighbourly polytope
bool isNeighbourly(vector<PointConfig>& c){
  bool result = true;
  for(int i = 0; i < c.size(); ++i){
    for(int j = i+1; j < c.size(); ++j){
      c[i].consider = false;
      c[j].consider = false;
      result = result and configIntersects(c);
      c[i].consider = true;
      c[j].consider = true;
    }
  }
  return result;
}
//Enumerates, recursively, all possible point configurations that we consider.
  void enumeratePointConfigurations(vector<PointConfig>& c, int index){
    if(index == c.size()){
      if(isNeighbourly(c)){
        printConfig(c);
      }
    }else{
      while(c[index].positive and index < c.size()) ++index;
      if(index == c.size()) enumeratePointConfigurations(c,index);
      else{
        c[index].index_point = 0;
        enumeratePointConfigurations(c,index+1);
        c[index].index_point = 1;
        enumeratePointConfigurations(c,index+1);
      }
    }   
    }
//Checks if we have already visited a combinatorially equivalent positive/negative configuration
bool isRepeatedSignConfig(vector<PointConfig>& c, vector < vector<bool> >& table){
  bool result = false;
  for(int i = 0; i < table.size(); ++i){
    result = result or check_dihedral(c, table[i]);
  }
  return result;
}
//Enumerates, recursively, all possible sign configurations with potential to be neighbourly, i.e, only those
//consisting of 3 or 4 negative points
void enumerateSignConfigurations(vector<PointConfig>& c, vector< vector<bool> >& table, int numNegs, int index){
    if(numNegs < 0){
      return;
    }else if(index == c.size()){
     if(numNegs != 0) return;
     //printConfig(c);
     if(!isRepeatedSignConfig(c, table)) {
       enumeratePointConfigurations(c,0);
       vector<bool> aux = vector<bool>(8);
       for(int i = 0; i < aux.size(); ++i){
         aux[i] = c[i].positive;
       }
       table.push_back(aux);
     } 
    }else{
      c[index].positive = false;
      enumerateSignConfigurations(c, table, numNegs - 1, index+1);
      c[index].positive = true;
      enumerateSignConfigurations(c, table, numNegs, index+1);
    }
  }
//General method that starts the enumeration
void enumerateConfigurations(vector<PointConfig>& c, vector < vector<bool> >& table){
    enumerateSignConfigurations(c, table, 3, 0);
    enumerateSignConfigurations(c, table, 4, 0);
  }

int main(int argc, char *argv[]) {
  const double pi = acos(-1);
  const double ctn = sin(pi/4);
  const double offset = 0.1;
  const double innerp = 1/2.0 - offset;
  const double innerpa = ctn - offset;
  vector<Point> outterConfig = vector<Point>(8);
  vector<Point> innerConfig = vector<Point>(8);
  outterConfig[0] = Point(0,1);
  outterConfig[1] = Point(ctn,ctn);
  outterConfig[2] = Point(1,0);
  outterConfig[3] = Point(ctn,-ctn);
  outterConfig[4] = Point(0,-1);
  outterConfig[5] = Point(-ctn,-ctn);
  outterConfig[6] = Point(-1,0);
  outterConfig[7] = Point(-ctn,ctn);

  innerConfig[0] = Point(0,innerpa);
  innerConfig[1] = Point(innerp,innerp);
  innerConfig[2] = Point(innerpa,0);
  innerConfig[3] = Point(innerp,-innerp);
  innerConfig[4] = Point(0,-innerpa);
  innerConfig[5] = Point(-innerp,-innerp);
  innerConfig[6] = Point(-innerpa,0);
  innerConfig[7] = Point(-innerp,innerp);
  vector<PointConfig> galeConfig = vector<PointConfig>(8);
  galeConfig[0] = PointConfig(outterConfig[0], innerConfig[0]);
  galeConfig[1] = PointConfig(outterConfig[1], innerConfig[1]);
  galeConfig[2] = PointConfig(outterConfig[2], innerConfig[2]);
  galeConfig[3] = PointConfig(outterConfig[3], innerConfig[3]);
  galeConfig[4] = PointConfig(outterConfig[4], innerConfig[4]);
  galeConfig[5] = PointConfig(outterConfig[5], innerConfig[5]);
  galeConfig[6] = PointConfig(outterConfig[6], innerConfig[6]);
  galeConfig[7] = PointConfig(outterConfig[7], innerConfig[7]);
  vector< vector < bool > > lookupTable = vector< vector < bool > >();
  enumerateConfigurations(galeConfig, lookupTable);
  return 0;
}
