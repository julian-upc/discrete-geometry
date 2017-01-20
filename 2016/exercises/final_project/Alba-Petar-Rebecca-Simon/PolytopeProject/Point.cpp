//
//  Point.cpp
//  PolytopeProject
//


#include "Point.hpp"


Point::Point(vector<mpq_class> coord,bool sign):sign(sign){
    for(int i = 0; i<coord.size(); i++){
        mpq_class rational = coord[i];
        rational.canonicalize();
        coordinates.push_back(rational);
    }
}

Point::Point(vector<double> coord, bool sign):sign(sign){
    for(double dob : coord){
        mpq_class q(dob);
        coordinates.push_back(q);
    }
}

bool Point::get_sign() const{
    return sign;
}

void Point::set_sign(bool s){
	sign = s;
}

vector<mpq_class> Point::get_coordinates() const {return coordinates;}


bool Point::one_sign() const{
    int tempsgn = sgn(coordinates[0]);
    for(int i=1; i<coordinates.size(); i++){
        if(tempsgn == sgn(coordinates[i])){
        } else return false;
    }
    return true;
}

Point point_to_vector(Point point){
    throw runtime_error("Still has to be implemented!!!");
}

Point Point::minus (const Point &b) const{
    vector<mpq_class> coord(coordinates.size());
    vector<mpq_class> bcoord = b.get_coordinates();
    for(int i=0; i<coordinates.size(); i++){
        coord[i] = coordinates[i]-bcoord[i];
        coord[i].canonicalize();
    }
    return Point(coord,true);
}

Point Point::multiply(const mpq_class &scalar) const{
    vector<mpq_class> nextPoints;
    for(int i = 0; i<coordinates.size(); i++){
        mpq_class nextCoo = coordinates[i]*scalar;
        nextCoo.canonicalize();
        nextPoints.push_back(nextCoo);
    }
    return Point(nextPoints,sign);
    
}

void Point::print() const{
    for(int i=0; i<coordinates.size(); i++){
        cout << coordinates[i].get_str() << " ";
    }
    cout << endl;
}
void Point::printLatex() const{
    for(int i=0; i<coordinates.size(); i++){
        mpq_class num =coordinates[i];
        
        cout << "$\\" << "frac{" <<num.get_num() << "}{" << num.get_den() << "}$" << " ";
    }
    cout << endl;
}

