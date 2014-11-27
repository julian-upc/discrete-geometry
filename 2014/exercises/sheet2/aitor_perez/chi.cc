#include <string>
#include <iostream>
#include <vector>
#include "io.cc"

using namespace std;

int factorial(int n) {
  if (n == 0) return 1;
  return n*factorial(n-1);
}

Mat permutations(int n) {
  Mat perm(factorial(n), Vec(n));
  
  if (n == 1) {
    perm[0][0] = 0;
    return perm;
  }
  
  Mat subperm = permutations(n-1);
  for (int i=0; i<n; i++) {
    for (int j=0; j<subperm.size(); j++) {
      for (int k=0; k<n; k++) {
        if (k==0) {
          perm[i*subperm.size()+j][0] = i;
        } else {
          perm[i*subperm.size()+j][k] = subperm[j][k-1] < i ? subperm[j][k-1] : subperm[j][k-1] + 1;
        }
      }
    }
  }
  return perm;
}

Mat permute(const Mat &m, const Vec &perm) {
  Mat p(m.size(), Vec(m[0].size()));
  
  for (int i=0; i<m.size(); i++) {
    for (int j=0; j<m[0].size(); j++) {
      p[i][j] = m[i][perm[j]];
    }
  }
  
  return p;
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cout << "--- Incorrect number of parameters. Usage: ./a.out number_of_test_case" << endl;
    return 0;
  }

  int n = 4;
  
  string pathP = "test/P" + (string)argv[1];
  string pathQ = "test/Q" + (string)argv[1];
  string pathSol = "test/solution" + (string)argv[1];
  string pathA = "test/A" + (string)argv[1];
  
  Mat P = transpose(read4x4(pathP.c_str()));
  Mat Q = transpose(read4x4(pathQ.c_str()));
  
  if (P.size() == 0 || Q.size() == 0) {
    cout << "--- Error reading input: Matrices are smaller than expected or number size overflows the capacity of the compiler" << endl;
    return 0;
  }
  
  long long int detP = det(P);
  long long int detQ = det(Q);
  
  // Check if determinants are zero
  if (detP == 0 || detQ == 0) {
    cout << "--- Vertices are not affinely independent" << endl;
    return 0;
  }
  
  // Check if determinants differ in sign
  if (detP != detQ && detP != -detQ) {
    cout << "--- There does not exist any transformation, since det(P) = " << detP << " and det(Q) = " << detQ << endl;
    return 0;
  }
  
  // Iterate over all permutations of the columns of P
  vector<Mat> mats;
  Mat perm = permutations(4);
  for (int i=0; i<perm.size(); i++) {
    Mat Pperm = permute(P, perm[i]);
    long long int detPperm = det(Pperm);
    
    Mat A = divide(Q*adjugate(Pperm), detPperm);
    
    if (A.size() > 0 && !belongs(A, mats)) {
      mats.push_back(A);
      write(A, pathSol.c_str());
    }
  }
  
  // If no matrix was written, write the zero matrix
  if (mats.size() == 0) {
    Mat z = zeros(4);
    mats.push_back(z);
    write(z, pathSol.c_str());
  }

  // Check if the test matrix is in the solution vector
  Mat ATest = read4x4(pathA.c_str());
  if (ATest.size() > 0) {
    if (belongs(ATest, mats)) {
      cout << "--- " << mats.size() << " transformation(s) found. Provided test matrix is in the set [TEST SUCCESSFUL]." << endl;
    } else {
      cout << "--- " << mats.size() << " transformation(s) found. Provided test matrix is NOT in the set [TEST FAILED]." << endl;
    }
  }
}
