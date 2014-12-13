#include <iostream>
#include <vector>

using namespace std;

typedef double Base;
typedef vector<Base> Vec;
typedef vector<Vec> Mat;


void show(const Vec& v) {
  for(int i=0; i<v.size(); i++) {
    if (i != 0) {
      cout << " ";
    }
    cout << v[i];
  }
  cout << endl;
}

void show(const Mat& m) {
  for(int i=0; i<m.size(); i++) {
    show(m[i]);
  }
}

Mat A(int n) {
  Mat m(n-1, Vec(n+1, 0));
  
  for (int i=0; i<n-1; i++) {
    m[i][i+1] = 1;
    m[i][i+2] = -1;
  }
  
  return m;
}

Mat D(int n) {
  Mat m(n, Vec(n+1, 0));
  
  for (int i=0; i<n; i++) {
    if (i < n-1) {
      m[i][i+1] = 1;
      m[i][i+2] = -1;
    } else {
      m[n-1][n-1] = 1;
      m[n-1][n] = 1;
    }
  }
  
  return m;
}

Mat E(int n) {
  Mat m(n, Vec(n+1, 0));
  
  for (int i=0; i<n; i++) {
    if (i < n-2) {
      m[i][i+1] = 1;
      m[i][i+2] = -1;
    } else if (i == n-2){
      m[n-2][n-2] = 1;
      m[n-2][n-1] = 1;
    } else {
      for (int j=1; j<n+1; j++) {
        m[n-1][j] = 1;
      }
    }
  }
  
  return m;
}

Mat F(int n) {
  Mat m(n, Vec(n+1, 0));
  
  for (int i=0; i<n; i++) {
    if (i == n-1) {
      for (int j=1; j<n+1; j++) {
        m[n-1][j] = 1;
      }
    } else {
      m[i][i+1] = 1;
      m[i][i+2] = (i == n-2) ? 0 : -1;
    }
  }
  
  return m;
}

const Base tau = 1.6180339887498948482;

Mat H(int n) {
  Mat m(n, Vec(n+1, 0));
  
  if (n == 3) {
    m[0][1] = tau + 1;
    m[0][2] = -1;
    m[0][3] = tau;
    
    m[1][1] = 1;
    m[2][3] = 1;
  } else if (n == 4) {
    m[0][1] = tau + 1;
    m[0][2] = -1;
    m[0][3] = -1;
    m[0][4] = -1;
    
    for (int i=1; i<n; i++) {
      m[i][i] = -1;
      m[i][i+1] = 1;
    }
  }
  
  return m;
}

Mat generators(char c, int n) {
  if (c == 'A') {
    return A(n+1);
  } else if (c == 'D') {
    return D(n);
  } else if (c == 'E') {
    return E(n);
  } else if (c == 'F') {
    return F(n);
  } else if (c == 'H') {
    return H(n);
  } else {
    return Mat(0);
  }
}

