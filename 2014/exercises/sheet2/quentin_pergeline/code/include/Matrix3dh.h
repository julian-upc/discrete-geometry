#ifndef MATRIX3DH_H
#define MATRIX3DH_H

#include <Rational.h>
#include <array>
#include <fstream>
#include <sstream>

namespace Custom{

class Matrix3dh
{
    public:

        //Constructors
        Matrix3dh();
        Matrix3dh(const Matrix3dh &m);
        Matrix3dh(std::ifstream &infile);

        //Destructor
        virtual ~Matrix3dh();

        Matrix3dh &operator= (const Matrix3dh &rhs);
        const Rational operator[] (int i)const;
        Rational &operator[] (int i);
        Matrix3dh operator*(Matrix3dh &other) const;

        Rational inv(Matrix3dh &invOut) const;

        Matrix3dh& swapCol (int i, int j);

        bool isInteger()const;

        bool operator== (const Matrix3dh &other) const;
        void permute(int i,std::array<Matrix3dh, 24> &T, int *n);
        std::array<Rational, 16> getM()const;

        static Matrix3dh nulMatrix();
        static Matrix3dh idMatrix();
    protected:
    private:
        std::array<Rational, 16> M;
};
std::ostream& operator<<(std::ostream& os, const Matrix3dh& m);

class MyHash
{
friend class Matrix3dh;
public:
    std::size_t make_hash(const Rational& v) const
    {
        return std::hash<int>()(v.numerator())^std::hash<int>()(v.denominator());
    }

    void hash_combine(std::size_t& h, const std::size_t& v) const
    {
        h ^= v + 0x9e3779b9 + (h << 6) + (h >> 2);
    }

    std::size_t operator()(const Matrix3dh &v) const
    {
        std::size_t h =0;
        for (const Rational &i : v.getM()){
            hash_combine(h, make_hash(i));
        }
        return h;
    }
};

}
#endif // MATRIX3DH_H
