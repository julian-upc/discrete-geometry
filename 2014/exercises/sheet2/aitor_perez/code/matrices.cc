#include <iostream>
#include <vector>

using namespace std;

typedef vector<long long int> Vec;
typedef vector<Vec> Mat;


Vec operator*(const Mat &m, const Vec &v) {
  int n1 = m.size();
  int n2 = v.size();
  
  Vec prod(n1);
  
  for (int i=0; i<n1; i++) {
    prod[i] = 0;
    for (int j=0; j<n2; j++) {
      prod[i] += m[i][j]*v[j];
    }
  }
  
  return prod;
}

Mat operator*(const Mat &m1, const Mat &m2) {
  int n1 = m1.size();
  int n2 = m2.size();
  if (n1 < 1 || n2 < 1) {
    return Mat(0);
  }
  
  int n3 = m2[0].size();
  if (n3 < 1 || m1[0].size() != n2) {
    return Mat(0);
  }

  Mat prod(n1, Vec(n3));
  
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n3; j++) {
      prod[i][j] = 0;
      for (int k=0; k<n2; k++) {
        prod[i][j] += m1[i][k]*m2[k][j];
      }
    }
  }
  
  return prod;
}

Mat identity(int n) {
  Mat m(n, Vec(n));
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      m[i][j] = i == j ? 1 : 0;
    }
  }
  
  return m;
}

Mat zeros(int n) {
  Mat m(n, Vec(n));
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      m[i][j] = 0;
    }
  }
  
  return m;
}

Mat transpose(const Mat &m) {
  int n1 = m.size();
  if (n1 < 1) {
    return Mat(0);
  }
  int n2 = m[0].size();
  Mat t(n2, Vec(n1));
  
  for (int i=0; i<n2; i++) {
    for (int j=0; j<n1; j++) {
      t[i][j] = m[j][i];
    }
  }
  
  return t;
}

// Returns the matrix m without row i and column j
Mat minor(const Mat &m, int i, int j) {
  int n1 = m.size();
  if (n1 < 1 || i>=n1) {
    return m;
  }
  int n2 = m[0].size();
  if (n2 < 1 || j>=n2) {
    return m;
  }
  
  Mat x(n1-1, Vec(n2-1));
  
  for (int k=0; k<n1-1; k++) {
    int k2 = k < i ? k : k+1;
    for (int l=0; l<n2-1; l++) {
      int l2 = l < j ? l : l+1;
      x[k][l] = m[k2][l2];
    }
  }
  
  return x;
}

long long int det(const Mat &m) {
  int n = m.size();
  if (n<1 || m[0].size() != n) {
    return 0;
  }
  
  if (n == 1) {
    return m[0][0];
  }
  
  if (n == 2) {
    return m[0][0]*m[1][1] - m[0][1]*m[1][0];
  }
  
  long long int d = 0;
  for (int i=0; i<n; i++) {
    d += m[i][0] == 0 ? 0 : (1-2*(i%2))*m[i][0]*det(minor(m, i, 0));
  }
  
  return d;
}

Mat adjugate(const Mat &m) {
  int n1 = m.size();
  if (n1 < 1) {
    return Mat(0);
  }
  int n2 = m[0].size();
  Mat a(n2, Vec(n1));
  
  for (int i=0; i<n1; i++) {
    for (int j=0; j<n2; j++) {
      a[i][j] = (1-2*((i+j)%2))*det(minor(m, i, j));
    }
  }
  
  return transpose(a);
}

Mat divide(const Mat &m, long long int d) {
  Mat div(m.size(), Vec(m[0].size()));
  
  for (int i=0; i<m.size(); i++) {
    for (int j=0; j<m[0].size(); j++) {
      if (m[i][j]%d == 0) {
        div[i][j] = m[i][j]/d;
      } else {
        return Mat(0);
      }
    }
  }
  
  return div;
}

bool equals(const Mat &m1, const Mat &m2) {
  if (m1.size() != m2.size() || m1[0].size() != m2[0].size()) {
    return false;
  }
  
  for (int i=0; i<m1.size(); i++) {
    for (int j=0; j<m1[0].size(); j++) {
      if (m1[i][j] != m2[i][j]) {
        return false;
      }
    }
  }
  
  return true;
}

bool belongs(const Mat &m, const vector<Mat> &v) {
  for (int i=0; i<v.size(); i++) {
    if (equals(m, v[i])) {
      return true;
    }
  }
  
  return false;
}

/*
Mat inverse(const Mat &m) {
  int n = m.size();
  if (n<1 || m[0].size() != n) {
    return Mat(0);
  }
  Mat inv(n, Vec(n));
  
  long long int d = det(m);
  Mat a = adjugate(m);
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (a[i][j]%d == 0) {
        inv[i][j] = a[i][j]/d;
      } else {
        return Mat(0);
      }
    }
  }
  
  return inv;
}
*/
