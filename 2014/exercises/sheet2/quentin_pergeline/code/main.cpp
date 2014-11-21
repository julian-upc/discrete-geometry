#include <unordered_set>
#include <iostream>
#include <Rational.h>
#include <array>
#include <Matrix3dh.h>
#include <set>
#include <fstream>

using Custom::Rational;
using Custom::Matrix3dh;
using Custom::MyHash;

namespace Custom{

void displayMatrix(Rational M[16]){
    for (int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            std::cout<<M[4*i+j].numerator()<<"/"<<M[4*i+j].denominator()<<" ";
        }
        std::cout<<std::endl;
    }
}

bool getMatrix(std::string name, Matrix3dh &P){
    std::ifstream file(name,std::ios::in);
    if (file){
        P = Matrix3dh(file);
        return true;
    } else {
        return false;
    }
}

void writeOutputFile(std::string name, std::unordered_set<Matrix3dh, MyHash> &Sol, Matrix3dh &P, Matrix3dh &Q){
    std::ofstream ofs (name, std::ofstream::trunc);
    if (ofs){
        ofs<<"============================= Input 1 ============================="<<std::endl;
        ofs<<P;
        ofs<<"============================= Input 2 ============================="<<std::endl;
        ofs<<Q;
        ofs<<"============================= Solution(s) ========================="<<std::endl;
        for(Matrix3dh s : Sol){
            ofs<< s <<std::endl;
        }
        ofs.close();
    } else {
        throw std::invalid_argument("Impossible to write the file");
    }
}

std::unordered_set<Matrix3dh, MyHash> Chi(Matrix3dh P, Matrix3dh Q){

    Matrix3dh Pinv;
    Rational detP = P.inv(Pinv);
    Rational detQ = Q.inv(Pinv);
    if (detP == 0 || detQ == 0)
        throw std::invalid_argument("Not 3-dimensional polytope");

    std::unordered_set<Matrix3dh, MyHash> Sol;
    std::array<Matrix3dh, 24> T;
    int *n = new int;
    *n = 0;
    P.permute(0,T,n);

    Matrix3dh sol;
    Rational det;
    for (int i=0; i<24; i++){
        sol = Q*T[i];
        if (sol.isInteger()){
            det = sol.inv(Pinv);
            if (det == 1 || det == -1){
                Sol.insert(sol);
            }
        }
    }
    if (Sol.empty()){
        Sol.insert(Matrix3dh().nulMatrix());
    }
    return Sol;
}

}

int main(int argc,  char** argv){
    std::string name_input1 = "test/input/matrix_P1.txt";
    std::string name_input2 = "test/input/matrix_P2.txt";
    std::string name_output = "test/output/matrix_P1_P2.txt";

    if (argc==4) {
        name_input1 = argv[1];
        name_input2 = argv[2];
        name_output = argv[3];
    } else if (argc != 1) {
        std::cerr<<"Wrong argument number: main name_input1 name_input2 name_output\n or main"<<std::endl;
        return -1;
    }

    Matrix3dh P,Q;
    std::unordered_set<Matrix3dh, MyHash> Sol;

    try {
        if (getMatrix(name_input1,P) && getMatrix(name_input2,Q)) {
            Sol = Chi(P,Q);
            writeOutputFile(name_output, Sol, P, Q);
        } else {
            std::cerr<<"Unable to read both files";
        }
    }
    catch(std::exception &e) {
        std::cerr << "Error: " << e.what() <<std::endl;
        return -1;
    }
    return 0;
}


