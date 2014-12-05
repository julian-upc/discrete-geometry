#include "Matrix3dh.h"
#include <stdexcept>
#include <limits>

namespace Custom{

Matrix3dh::Matrix3dh()
{
    M = {0};
}

Matrix3dh::Matrix3dh(const Matrix3dh &m){
    *this = m;
}

Matrix3dh::Matrix3dh(std::ifstream &infile){
   std::string line;
    int j(0);
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        long long int a;
        for(int i = 0; i<4; i++){
            if (!(iss >> a) || iss.fail()) {throw std::invalid_argument("Invalid input");}

            if (i==0 && a !=1) throw std::invalid_argument("Non homogeneous coordinates");
            else {
                M[4*i+j] = a;
            }
        }

        if (iss >> a) { throw std::invalid_argument("Too much data on a line");}
        j++;

        if (j>5)//Too much lines!
            throw std::invalid_argument("Too much lines");
    }
    if (j != 4)
        throw std::invalid_argument("Not enough lines");
}


Matrix3dh::~Matrix3dh()
{
    //dtor
}

std::array<Rational, 16> Matrix3dh::getM() const {return M;}

Matrix3dh &Matrix3dh::operator= (const Matrix3dh &rhs){
     for (int i=0; i<4;i++)
        for (int j=0; j<4;j++)
            (*this)[4*i+j] = rhs[4*i+j];

    return *this;
}

 bool Matrix3dh::operator== (const Matrix3dh &other) const{
    bool eq = true;
    for (int i=0; i<4;i++)
        for (int j=0; j<4;j++){
            eq = eq&&((*this)[i*4+j]==other[4*i+j]);
        }

    return eq;
}

const Rational Matrix3dh::operator[] (int i)const{
    if (i>15 || i<0) throw 0;
    const Rational r = this->M[i];
    return r;
}

Rational& Matrix3dh::operator[] (int i){
    if (i>15 || i<0) throw 0;
    return this->M[i];
}

std::ostream& operator<<(std::ostream& os, const Matrix3dh& m){
    for (int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            os<<m[4*i+j];
        }
        os<<std::endl;
    }
    return os;
}
 Matrix3dh Matrix3dh::operator*(Matrix3dh &other) const{
    Matrix3dh prod;
    Rational t = 0;

    for (int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            for(int k=0; k<4; k++){
                t = t + (*this)[4*i+k]*other[4*k+j];
            }
            prod[4*i+j] = t;
            t = 0;
        }
    }
    return prod;
}

bool Matrix3dh::isInteger()const{
    bool eq = true;
    for (int i=0; i<4;i++)
        for (int j=0; j<4;j++){
            eq = eq&&((*this)[i*4+j].isInteger());
        }

    return eq;
}
/*  
 * Compute determinant. If non-zero, return its value and the parameter henceworth stocks the inverse matrix. 
 */
Rational Matrix3dh::inv(Matrix3dh &invOut) const{
    Rational inv[16], det;
    std::array<Rational, 16> m = this->M;

    inv[0] = m[5]  * m[10] * m[15] -
             m[5]  * m[11] * m[14] -
             m[9]  * m[6]  * m[15] +
             m[9]  * m[7]  * m[14] +
             m[13] * m[6]  * m[11] -
             m[13] * m[7]  * m[10];

    inv[4] = -m[4]  * m[10] * m[15] +
              m[4]  * m[11] * m[14] +
              m[8]  * m[6]  * m[15] -
              m[8]  * m[7]  * m[14] -
              m[12] * m[6]  * m[11] +
              m[12] * m[7]  * m[10];

    inv[8] = m[4]  * m[9] * m[15] -
             m[4]  * m[11] * m[13] -
             m[8]  * m[5] * m[15] +
             m[8]  * m[7] * m[13] +
             m[12] * m[5] * m[11] -
             m[12] * m[7] * m[9];

    inv[12] = -m[4]  * m[9] * m[14] +
               m[4]  * m[10] * m[13] +
               m[8]  * m[5] * m[14] -
               m[8]  * m[6] * m[13] -
               m[12] * m[5] * m[10] +
               m[12] * m[6] * m[9];

    inv[1] = -m[1]  * m[10] * m[15] +
              m[1]  * m[11] * m[14] +
              m[9]  * m[2] * m[15] -
              m[9]  * m[3] * m[14] -
              m[13] * m[2] * m[11] +
              m[13] * m[3] * m[10];

    inv[5] = m[0]  * m[10] * m[15] -
             m[0]  * m[11] * m[14] -
             m[8]  * m[2] * m[15] +
             m[8]  * m[3] * m[14] +
             m[12] * m[2] * m[11] -
             m[12] * m[3] * m[10];

    inv[9] = -m[0]  * m[9] * m[15] +
              m[0]  * m[11] * m[13] +
              m[8]  * m[1] * m[15] -
              m[8]  * m[3] * m[13] -
              m[12] * m[1] * m[11] +
              m[12] * m[3] * m[9];

    inv[13] = m[0]  * m[9] * m[14] -
              m[0]  * m[10] * m[13] -
              m[8]  * m[1] * m[14] +
              m[8]  * m[2] * m[13] +
              m[12] * m[1] * m[10] -
              m[12] * m[2] * m[9];

    inv[2] = m[1]  * m[6] * m[15] -
             m[1]  * m[7] * m[14] -
             m[5]  * m[2] * m[15] +
             m[5]  * m[3] * m[14] +
             m[13] * m[2] * m[7] -
             m[13] * m[3] * m[6];

    inv[6] = -m[0]  * m[6] * m[15] +
              m[0]  * m[7] * m[14] +
              m[4]  * m[2] * m[15] -
              m[4]  * m[3] * m[14] -
              m[12] * m[2] * m[7] +
              m[12] * m[3] * m[6];

    inv[10] = m[0]  * m[5] * m[15] -
              m[0]  * m[7] * m[13] -
              m[4]  * m[1] * m[15] +
              m[4]  * m[3] * m[13] +
              m[12] * m[1] * m[7] -
              m[12] * m[3] * m[5];

    inv[14] = -m[0]  * m[5] * m[14] +
               m[0]  * m[6] * m[13] +
               m[4]  * m[1] * m[14] -
               m[4]  * m[2] * m[13] -
               m[12] * m[1] * m[6] +
               m[12] * m[2] * m[5];

    inv[3] = -m[1] * m[6] * m[11] +
              m[1] * m[7] * m[10] +
              m[5] * m[2] * m[11] -
              m[5] * m[3] * m[10] -
              m[9] * m[2] * m[7] +
              m[9] * m[3] * m[6];

    inv[7] = m[0] * m[6] * m[11] -
             m[0] * m[7] * m[10] -
             m[4] * m[2] * m[11] +
             m[4] * m[3] * m[10] +
             m[8] * m[2] * m[7] -
             m[8] * m[3] * m[6];

    inv[11] = -m[0] * m[5] * m[11] +
               m[0] * m[7] * m[9] +
               m[4] * m[1] * m[11] -
               m[4] * m[3] * m[9] -
               m[8] * m[1] * m[7] +
               m[8] * m[3] * m[5];

    inv[15] = m[0] * m[5] * m[10] -
              m[0] * m[6] * m[9] -
              m[4] * m[1] * m[10] +
              m[4] * m[2] * m[9] +
              m[8] * m[1] * m[6] -
              m[8] * m[2] * m[5];

    det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];

    if (det == 0)
        return det;

    det = 1 / det;
    int i;
    for (i = 0; i < 16; i++)
        invOut[i] = inv[i] * det;

    return det;
}

Matrix3dh& Matrix3dh::swapCol (int i, int j){
    if ((i<0 || i>3) && (j<0 || j>3)) throw 0;
    if (i==j) return *this;

    Rational temp[4];
    for(int k = 0; k<4; k++){
        temp[k] = (*this)[4*k+i];
        (*this)[4*k+i] = (*this)[4*k+j];
    }
    for(int k = 0; k<4; k++){
        (*this)[4*k+j] =  temp[k];
    }
    return *this;
}

Matrix3dh Matrix3dh::nulMatrix(){
    Matrix3dh nul = Matrix3dh();
    nul[0]=0;
    nul[5]=0;
    nul[10]=0;
    nul[15]=0;
    return nul;
}

Matrix3dh Matrix3dh::idMatrix(){
    Matrix3dh id = Matrix3dh();
    id[0]=1;
    id[5]=1;
    id[10]=1;
    id[15]=1;
    return id;
}

/*
* Computes the 24 prmutations of the columns and stock then in T
*/
void Matrix3dh::permute(int i,std::array<Matrix3dh, 24> &T, int *n){
    int j;
    Matrix3dh m;
    if (i == 4){
        this->inv(m);
        T[*n] = m;
        (*n)++;
    }
    else{
        for (j = i; j < 4; j++)
        {
            swapCol(i, j);
            permute(i+1,T,n);
            swapCol(i, j);

        }
    }
}

}
