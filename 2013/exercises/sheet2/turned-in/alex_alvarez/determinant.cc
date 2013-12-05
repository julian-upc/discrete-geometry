#include <iostream>
#include <vector>
#include <gmpxx.h>
#include <cstdlib>
#include <set>
#include <map>
#include <algorithm>
using namespace std;

typedef long long ll;
typedef mpq_class T;

vector<T> operator*(const T& x, const vector<T>& v) {
  vector<T> res(v.size());
  for (int i = 0; i < v.size(); ++i) res[i] = x*v[i];
  return res;
}

vector<T> operator-(const vector<T>& x, const vector<T>& y) {
  vector<T> res(x.size());
  for (int i = 0; i < x.size(); ++i) res[i] = x[i] - y[i];
  return res;
}

void operator/=(vector<T>& v, const T& x) {
  //Done in that order to be able to pass the x by reference
  for (int i = v.size() - 1; i >= 0; --i) v[i] /= x;
}

T compute_determinant(vector<vector<T> >& mat) {
  int r = mat.size(), c = mat[0].size();
  T det(1, 1);

  for (int i = 0; i < r; ++i)
    for (int j = 0; j < c; ++j)
      mat[i][j].canonicalize();

  for (int i = 0, j = 0; i < r and j < c; ++j) {
    int nz = i;
    while (nz < r and mat[nz][j] == 0) ++nz;
    if (nz == r) return T(0);
    if (nz != i) det = -det;
    swap(mat[i], mat[nz]);
    det *= mat[i][j];
    mat[i] /= mat[i][j];
    for (int k = i + 1; k < r; ++k)
      mat[k] = mat[k] - mat[k][j]*mat[i];
    ++i;
  }
  return det;
}

int main(int argc, char** argv) {
  
  if (argc < 4) {
    cout << "usage: program limit max_dim seed" << endl;
    cout << "\t\tlimit: number of random matrices to test for each dimension" << endl;
    cout << "\t\tmax_dim: maximum dimension to perform the test" << endl;
    cout << "\t\tseed: a seed for the randon number generator" << endl;
    return 1;
  }
  
  const int LIMIT = atoi(argv[1]);
  const int MAX_DIM = atoi(argv[2]);
  srand(atoi(argv[3]));
  
  for (int d = 3; d <= MAX_DIM; ++d) {
    
    map<T, int> hist;
    vector<vector<T> > mat(d + 1, vector<T>(d + 1, 1));
        
    for (int step = 0; step < LIMIT; ++step) {      
      
      for (int i = 0; i < d + 1; ++i) {
        
        mat[i][0] = 1;
        for (int j = 0; j < d; ++j) {
          mat[i][j + 1] = (rand() & 1 ? 1 : -1);
        }
      }
      
      ++hist[abs(compute_determinant(mat))];
    }
    
    cout << "Results for dimension " << d << ":" << endl;
    for (map<T, int>::iterator it = hist.begin(); it != hist.end(); ++it) {
      cout << it->first << " " << it->second << endl;    
    }

  }
}
