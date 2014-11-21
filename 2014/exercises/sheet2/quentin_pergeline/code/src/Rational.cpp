#include "Rational.h"
#include <iostream>

namespace Custom{

Rational::Rational()
{
    n=0;
    d=1;
}

Rational::Rational(long long int n){
    this->n = n;
    this->d = 1;
}

/*
 * Called whenever an operation on Rational is made. It prevents integer overflow.
*/
void Rational::optimize(){
    long long int n1 = this->n;
    long long int d1 = this->d;

    if (d1==0) throw 0;
    else if(d1<0){
        n1 = -n1;
        d1 = -d1;
    }

    if (n1%d1==0){
        n1 = n1/d1;
        d1 = 1;
    } else {
        long long int g = gcd(n1,d1);
        if (g>1){
            n1 = n1/g;
            d1 = d1/g;
        }
    }
    this->n = n1;
    this->d = d1;
}

Rational::Rational(long long int n, long long int d)
{

    this->n = n;
    this->d = d;

    this->optimize();
}

Rational::Rational (const Rational & r){
    this->n = r.n;
    this->d = r.d;
}

Rational::~Rational()
{
    //dtor
}

/*
 *Input: integer
 *Output: absolute value of the integer
 */
long long int Rational::abs(long long int a){
    if (a<0) return -a;
    return a;
}

/*
 *Input: two integers a and b
 *Output: the positive representant of the greater common divisor of a and b.
 */
long long int Rational::gcd(long long int a, long long int b){
    a = abs(a);
    b = abs(b);
    long long int t = 0;
    while (b != 0){
        t = b;
        b = a%b;
        a = t;
    }
    return a;
}

long long int Rational::numerator()const{
    return this->n;
}
long long int Rational::denominator()const{
    return this->d;
}

Rational & Rational::operator= (const Rational &rhs){
    this->n = rhs.n;
    this->d = rhs.d;

    return *this;
}

const Rational  Rational::operator+ (const Rational &other)const{
    Rational sum(this->n * other.d + other.n * this->d, this->d * other.d);

    sum.optimize();

    return sum;
}


const Rational  Rational::operator-()const {
    Rational opp(-this->n, this->d);
    return opp;
}

const Rational  Rational::operator- (const Rational &other)const{

    return *this + -other;
}

const Rational Rational::operator*(const Rational &other)const{
    Rational prod(this->n * other.n , this->d * other.d);
    prod.optimize();
    return prod;
}

const Rational Rational::operator/(const Rational &other)const{
    Rational quot(this->n * other.d, this->d * other.n);

    quot.optimize();
    return quot;
}

bool Rational::operator== (const Rational &other)const{
    return (this->n * other.d == other.n * this->d);
}
bool Rational::operator!= (const Rational &other) const{
    return !(*this==other);
}

bool Rational::isInteger() const{
    return (denominator()==1);
}

std::ostream& operator<<(std::ostream& os, const Rational& r)
{
    if (r.denominator()==1)
        os <<" "<< r.numerator() <<" ";
    else
        os << r.numerator() << "/"<< r.denominator()<<" ";
    return os;
}

Rational operator+ (const long long int i,  Rational r){
    return Rational(i)+r;
}

Rational operator- (const  long long int i,  Rational r){
    return Rational(i)-r;
}
Rational operator* (const  long long int i,  Rational r){
    return Rational(i)*r;
}
Rational operator/ (const  long long int i,  Rational r){
    return Rational(i)/r;
}

}
