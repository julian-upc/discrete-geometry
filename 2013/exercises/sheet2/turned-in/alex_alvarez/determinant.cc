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

T inverse(const T& x) {
  T res(x.get_den(), x.get_num());
  res.canonicalize();
  return res;
}

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

T gauss(vector<vector<T> >& mat) {
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
    det /= inverse(mat[i][j]);
    mat[i] = inverse(mat[i][j])*mat[i];
    for (int k = i + 1; k < r; ++k)
      mat[k] = mat[k] - mat[k][j]*mat[i];
    ++i;
  }
  return det;
}

ll random_int(int d) {
  ll x = 0;
  for (int i = 0; i < d; ++i) if (rand()&1) x |= (1 << i);
  return x;
}

int main() {
  srand(time(NULL));
  const int LIMIT = 1000000;
  const int MAX_DIM = 10;
  
  for (int d = 3; d <= MAX_DIM; ++d) {
    
    map<T, int> hist;
        
    for (int step = 0; step < LIMIT; ++step) {
      set<ll> points;
      while (points.size() < d + 1) points.insert(random_int(d));
      
      vector<vector<T> > mat(d + 1, vector<T>(d + 1, 1));
      
      set<ll>::iterator it = points.begin();
      for (int i = 0; i < d + 1; ++i) {
        ll x = *it;
        
        for (int j = 0; j < d; ++j){
          if ((x >> j)&1) mat[i][j + 1] = 1;
          else mat[i][j + 1] = -1;
        }
        
        ++it;
      }
      T aux = gauss(mat);
      if (aux < 0) aux = -aux;
      ++hist[aux];
    }
    
    cout << "Results for dimension " << d << ":" << endl;
    for (map<T, int>::iterator it = hist.begin(); it != hist.end(); ++it) {
      cout << it->first << " " << it->second << endl;    
    }

  }
}
