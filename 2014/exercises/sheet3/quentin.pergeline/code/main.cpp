#include <iostream>
#include <unordered_set>
#include "symbolicc++.h"

//The unordered_set demanded me to define a custom hash in order to store the orbits.
class MyHash
{
friend class Vector<Symbolic>;
public:
    std::size_t make_hash(const Symbolic& v) const
    {
        return std::hash<double>()(v);
    }

    void hash_combine(std::size_t& h, const std::size_t& v) const
    {
        h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
    }

    std::size_t operator()(const Vector<Symbolic> &v) const
    {
        std::size_t h =0;
        for (const Symbolic &i : v){
            hash_combine(h, make_hash(i));
        }
        return h;
    }
};

// Any point is a Pt_t
typedef Vector<Symbolic> Pt_t;

//SymbolicC++ did not implement the scalar product lambda * v. Where lambda is a real and v is a point.
Pt_t operator* (const Symbolic k,  const Pt_t &r);

// I am forced to use doubles to check unicity of the orbits. This check uses a distance between points.
double dist(Pt_t p, Pt_t  q);

// Compute the reflection of the point p by the hyperphlane whose one normal vector is norm
inline Pt_t reflection(Pt_t  norm, Pt_t  p){return p - 2* ( (norm | p) / (norm | norm) )*norm;}

// Computes all the orbits of initial point p,  by the reflection group whose generators are normList.
void enumOrbits(Pt_t &p, std::size_t prev, vector<Pt_t> &normList, std::unordered_set<Vector<Symbolic>, MyHash> &enumerated, int *n);

// Compute an initial point for the orbits, in case the user doesn't give one.
Pt_t initial_orbit (char X, int n);

// The final function called in the main. first is by default an initial point not lying in any hyperplane.
int orbit(char X, int n, Pt_t first);



// The initial point used, depending on the reflection group used, was made up "by hand".
Pt_t initial_orbit (char X, int n){
    Pt_t v(n+1);
    for (int i =0; i < n+1; i++){
        v[i] = i+1;
    }

    if ((X == 'F') ) {v[4] += 0.1;}
    else if (X == 'E') {v[2] += 0.1;}

    return v;
}

double dist(Pt_t p, Pt_t  q){
    if (p.size() != q.size()) throw std::invalid_argument("Vectors of different size.");
    double d = 0.0;
    double d1 = 0.0;
    double d2 = 0.0;
    for (size_t i =0; i < p.size(); i++){
        d1 = p[i];
        d2 = q[i];
        d += (d1 - d2)*(d1 - d2);
    }
    return d;
}

Pt_t operator* (const Symbolic k,  const Pt_t &r){
    Pt_t v;
    for (std::size_t i=0; i< r.size(); i++){
        v.push_back(k * r[i]);
    }
    return v;
}

const vector<Pt_t> generate_An (int n){

    if (n < 2) throw std::invalid_argument("An only implemented for n greater than 1");

    vector<Pt_t> normList;
    Pt_t v(n+1);

    for (int i = 0; i < n-1; i++){
        v[i+1] = 1;
        v[i+2] = -1;
        normList.push_back(v);
        v.reset(n+1);
    }
    return normList;
}

const vector<Pt_t> generate_Dn (int n){

    if (n < 2) throw std::invalid_argument("En only implemented for n greater than 1");

    vector<Pt_t> normList;
    Pt_t v(n+1);

    for (int i = 0; i < n-1; i++){
        v[i+1] = 1;
        v[i+2] = -1;
        normList.push_back(v);
        v.reset(n+1);
    }

    v[n-1] = 1;
    v[n] = 1;
    normList.push_back(v);

    return normList;
}

const vector<Pt_t> generate_En (int n){

    if (n < 2) throw std::invalid_argument("En only implemented for n greater than 1");

    vector<Pt_t> normList;
    Pt_t v(n+1);

    for (int i = 0; i < n -2; i++){
        v[i+1] = 1;
        v[i+2] = -1;
        normList.push_back(v);
        v.reset(n+1);
    }

    v[n-2] = 1;
    v[n-1] = 1;
    normList.push_back(v);

    Pt_t w(n+1, sqrt(2)/2);
    w[0] = 0;
    w = (1/sqrt(2))*w;
    normList.push_back(w);

    return normList;
}

const vector<Pt_t> generate_F4 (){
    vector<Pt_t> normList;

    Pt_t v1(5, 0), v2(5,0), v3(5,0), v4(5, -sqrt(2)/2 );
    v1[1] = 1;
    v1[2] = -1;

    v2[2] = 1;
    v2[3] = -1;

    v3[3] = 1;

    v4[0] = 0;
    v4 = (1/sqrt(2))*v4;

    normList.push_back(v1);
    normList.push_back(v2);
    normList.push_back(v3);
    normList.push_back(v4);

    return normList;
}

const vector<Pt_t> generate_Hn (int n){
    Symbolic tau =   (sqrt(5) + 1)/2;
    vector<Pt_t> normList;

    if (n == 3){
        Pt_t v1(4, 0), v3(4,0);
        v1[1]=2;
        v3[3]=2;
        Pt_t v2 (4);
        v2[1] = -tau;
        v2[2] = 1/tau;
        v2[3] = -1;

        normList.push_back(v1);
        normList.push_back(v3);
        normList.push_back(v2);
    } else if (n == 4){
        Symbolic tau =   (sqrt(5) + 1)/2;
        Pt_t v(5);
        v[0] = 4;
        v[1] = -tau;
        v[2] = 1/tau;
        v[3] = 1/tau;
        v[4] = 1/tau;
        normList.push_back(v);
        v.reset(5);
        for (int i = 1; i < 4; i++){
            v[i+1] = 1;
            v[i] = -1;
            normList.push_back(v);
            v.reset(5);
        }
    } else {
        throw std::invalid_argument("Only H3 and H4 are implemented.");
    }
    return normList;
}

/* THIS IS the recursive function used to compute the orbit.
 * For any point p, we compute the images of p by all the generators, except the one (prev) whose p is the image by a previous point.
 * Whenever an image is not already in the unordered set of discovered points called enumerated, we add it and call the recursive function 
 * taking the image point just added as the new p.
 **/
void enumOrbits(Pt_t &p, std::size_t prev, vector<Pt_t> &normList, std::unordered_set<Vector<Symbolic>, MyHash> &enumerated, int *n){
    Pt_t newPt;
    for (std::size_t i = 0; i< normList.size(); i++){
        if (i == prev){
             // Loss of time to compute back the points where we come from.
            continue;
        }
        newPt = reflection(normList[i], p);

        bool newOrbit = true;
        for (auto v : enumerated){
            if (dist(newPt,v) < (double)(0.0001) ) newOrbit = false;
        }
        if (newOrbit && *n < 30000){
            enumerated.insert(newPt);
            (*n) += 1;
            enumOrbits(newPt, i, normList, enumerated, n);
        }
    }
}

/* This is the main function called by program. It deals with setting the system with the right reflection group X, its index n, and a first point.
 * It then make the first call to the recursive function.
 */
int orbit(char X, int n, Pt_t first= Vector<Symbolic>()){
    vector<Pt_t> normList;

    switch (X) {
    case 'A':
        if (n>1){
            normList = generate_An(n);
        } else {
            throw std::invalid_argument("An only implemented for n greater than 1");
        }
        break;

    case 'H':
        if (n == 3 || n == 4){
            normList = generate_Hn(n);
        } else {
            throw std::invalid_argument("Only H3 and H4 are implemented.");
        }
        break;

    case 'F':
        if (n == 4 ){
            normList = generate_F4();
        } else {
            throw std::invalid_argument("Only F4 is implemented.");
        }
        break;

    case 'E':
        if (n > 1 ){
            normList = generate_En(n);
        } else {
            throw std::invalid_argument("En only implemented for n greater than 1");
        }
        break;

    case 'D':
        if (n > 1 ){
            normList = generate_Dn(n);
        } else {
            throw std::invalid_argument("En only implemented for n greater than 1");
        }
        break;

    default:
        throw std::invalid_argument("first argument must be either one of A, D, E, P, H.");
        break;
    }

    int *order =new int(1);
    std::unordered_set<Vector<Symbolic>, MyHash> Already;

    Pt_t v;
    if (first.size() == 0){
        v = initial_orbit(X, n);
    } else {
        v = first;
    }
    //for(Pt_t w : normList){v += w;}

    Already.insert(v);
    enumOrbits(v, normList.size(), normList, Already, order);

    return *order;
}

/*
 * I did not know how to let the user enter a personal first vector. So the only way is changing the code and recompiling (instructions to recompile in the diary)
 * But Symbolicc++ is far less easy to handle to put vectors than numpy in Python...
 */
int main(int argc, char *argv[]){


    char group_input = 'H';
    int order_input = '3'-'0';

    if (argc==3) {
        group_input = argv[1][0];
        order_input = argv[2][0] - '0';

    } else if (argc != 1) {
        std::cerr<<"Wrong argument(s) : orbit A|D|E|F|H  order\n"<<std::endl;
        return -1;
    }

    try {
        std::cout << "Group to analyse: " << group_input << std::endl;
        std::cout << "Order: " << order_input << std::endl;
        std::cout << "Computed number of orbits: " << orbit(group_input, order_input) << std::endl;
    }
    catch(std::exception &e) {
        std::cerr << "Error: " << e.what() <<std::endl;
        return -1;
    }
    return 0;
}

