#ifndef RATIONAL_H
#define RATIONAL_H

#include <iostream>

namespace Custom{

class Rational
{
    public:

        //Constructors
        Rational(); //numerator=0 and denominator=1
        Rational(long long int n); //numerator=n and denominator=1
        Rational(long long int n, long long int d);//normalized rational equivalent to n/d
        Rational (const Rational & r); // Copy constructor

        //Accessors
        long long int numerator() const;
        long long int denominator() const;

        //Destructor
        virtual ~Rational();

        //Operator Overloadings
        Rational &operator= (const Rational &rhs);
        const Rational operator+ (const Rational &other) const;
        const Rational operator-() const;
        const Rational operator- (const Rational &other) const;
        const Rational operator*(const Rational &other)const;
        const Rational operator/(const Rational &other) const;

        bool operator== (const Rational &other) const;
        bool operator!= (const Rational &other) const;

        bool isInteger() const;
    private:
        long long int n=0;
        long long int d=1;

        void optimize();
        static long long int gcd(long long int a, long long int b);
        static long long int abs(long long int a);
};

std::ostream& operator<<(std::ostream& os, const Rational& r);

Rational operator+ (const long long int i,  Rational r);
Rational operator* (const long long int i,  Rational r);
Rational operator/ (const long long int i,  Rational r);
Rational operator- (const long long int i,  Rational r);

}

#endif // RATIONAL_H
