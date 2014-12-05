#include <iostream>
#include <fstream>
#include "matrices.cc"

using namespace std;

Mat read4x4(const char* path) {
  int n = 4;
  Mat m(n, Vec(n));
  
  std::ifstream file(path);
  
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      if (!(file >> m[i][j])) {
        return Mat(0);
      }
    }
  }
  
  return m;
}

void write(const Mat &m, const char* path) {
  std::ofstream file(path, ios::app);
  
  for (int i=0; i<m.size(); i++) {
    for (int j=0; j<m[0].size(); j++) {
      if (j != 0) {
        file << " ";
      }
      file << m[i][j];
    }
    file << endl;
  }
}

void write(const vector<Mat> &v, const char* path) {
  for (int i=0; i<v.size(); i++) {
    write(v[0], path);
  }
}

void show(const Mat &m) {
  int n1 = m.size();
  if (n1 > 0) {
    int n2 = m[0].size();
    for (int i=0; i<n1; i++) {
      for (int j=0; j<n2; j++) {
        if (j != 0) {
          cout << " ";
        }
        cout << m[i][j];
      }
      cout << endl;
    }
  }
  cout << endl;
}
