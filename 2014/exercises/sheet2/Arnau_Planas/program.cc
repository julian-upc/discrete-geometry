#include<iostream>
#include<vector>
#include<fstream>
#include<string>

#define vi vector<long long>
#define mi vector<vi>

#define vd vector<double>
#define md vector<vd>

using namespace std;

mi P(24, vi(4));
vector<mi> list(24, mi(4, vi(4)));
vector<bool> used(24,0);

ofstream fs("output.txt");
ifstream fe("input.txt");


void initialize_P(){ //initalize a matrix with the permutations of four elements
P[0][0] = 0; P[0][1] = 1; P[0][2] = 2; P[0][3] = 3; //1
P[1][0] = 0; P[1][1] = 1; P[1][2] = 3; P[1][3] = 2; //-1
P[2][0] = 0; P[2][1] = 2; P[2][2] = 1; P[2][3] = 3; //-1
P[3][0] = 0; P[3][1] = 2; P[3][2] = 3; P[3][3] = 1; //1
P[4][0] = 0; P[4][1] = 3; P[4][2] = 1; P[4][3] = 2; //1
P[5][0] = 0; P[5][1] = 3; P[5][2] = 2; P[5][3] = 1; //-1
P[6][0] = 1; P[6][1] = 0; P[6][2] = 2; P[6][3] = 3; //-1
P[7][0] = 1; P[7][1] = 0; P[7][2] = 3; P[7][3] = 2; //1
P[8][0] = 1; P[8][1] = 2; P[8][2] = 0; P[8][3] = 3; //1
P[9][0] = 1; P[9][1] = 2; P[9][2] = 3; P[9][3] = 0; //-1
P[10][0] = 1; P[10][1] = 3; P[10][2] = 0; P[10][3] = 2; //-1 
P[11][0] = 1; P[11][1] = 3; P[11][2] = 2; P[11][3] = 0; //1
P[12][0] = 2; P[12][1] = 0; P[12][2] = 1; P[12][3] = 3; //1
P[13][0] = 2; P[13][1] = 0; P[13][2] = 3; P[13][3] = 1; //-1
P[14][0] = 2; P[14][1] = 1; P[14][2] = 0; P[14][3] = 3; //-1
P[15][0] = 2; P[15][1] = 1; P[15][2] = 3; P[15][3] = 0; //1
P[16][0] = 2; P[16][1] = 3; P[16][2] = 0; P[16][3] = 1; //1
P[17][0] = 2; P[17][1] = 3; P[17][2] = 1; P[17][3] = 0; //-1
P[18][0] = 3; P[18][1] = 0; P[18][2] = 1; P[18][3] = 2; //-1
P[19][0] = 3; P[19][1] = 0; P[19][2] = 2; P[19][3] = 1; //1
P[20][0] = 3; P[20][1] = 1; P[20][2] = 0; P[20][3] = 2; //1
P[21][0] = 3; P[21][1] = 1; P[21][2] = 2; P[21][3] = 0; //-1
P[22][0] = 3; P[22][1] = 2; P[22][2] = 0; P[22][3] = 1; //-1
P[23][0] = 3; P[23][1] = 2; P[23][2] = 1; P[23][3] = 0; //1
}

void read_matrix(mi & M, string s){ //reads a matrix
  const char *y = s.c_str();
  ifstream read(y);
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      read >> M[j][i];
    }
  }
  read.close();
}

void write_matrix(const mi & M){ //write a matrix
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      fs << M[i][j] << " ";
    }
    fs << endl;
  }
}

long long determinant(const mi & M){ //calculate the determinant of a matrix
  long long det1, det2, det3, det4, det;
  det1 = M[1][1]*M[2][2]*M[3][3] + M[2][1]*M[3][2]*M[1][3] + M[1][2]*M[2][3]*M[3][1] 
      - M[3][1]*M[2][2]*M[1][3] - M[3][2]*M[2][3]*M[1][1] - M[2][1]*M[1][2]*M[3][3];
      
  det2 = M[1][0]*M[2][2]*M[3][3] + M[2][0]*M[3][2]*M[1][3] + M[1][2]*M[2][3]*M[3][0] 
      - M[3][0]*M[2][2]*M[1][3] - M[3][2]*M[2][3]*M[1][0] - M[2][0]*M[1][2]*M[3][3];
      
  det3 = M[1][0]*M[2][1]*M[3][3] + M[2][0]*M[3][1]*M[1][3] + M[1][1]*M[2][3]*M[3][0] 
      - M[3][0]*M[2][1]*M[1][3] - M[3][1]*M[2][3]*M[1][0] - M[2][0]*M[1][1]*M[3][3];
  
  det4 = M[1][0]*M[2][1]*M[3][2] + M[2][0]*M[3][1]*M[1][2] + M[1][1]*M[2][2]*M[3][0] 
      - M[3][0]*M[2][1]*M[1][2] - M[3][1]*M[2][2]*M[1][0] - M[2][0]*M[1][1]*M[3][2];
      
  det = M[0][0]*det1 - M[0][1]*det2 + M[0][2]*det3 - M[0][3]*det4;
  return det;
}

mi minor(const mi & M, int k, int l){ //return the same matrix with zeros at row k, column l and 1 at position kl
  mi N(4, vi(4));
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      if (i == k or j == l) N[i][j] = 0;
      else N[i][j] = M[i][j];
      if (i == k and j == l) N[i][j] = 1;
    }
  }
  return N;
}

mi almost_inverse_matrix(const mi & M){ //return the inverse matrix multiplied by the determinant
  mi N(4, vi(4));
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      N[j][i] = determinant(minor(M, i, j));
    }
  }
  return N;
}

mi multiply(const mi & M, const mi & N){ //multiplies two matrices
  mi L(4, vi(4));
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      L[i][j] = 0;
      for (int k = 0; k < 4; k++) L[i][j] += M[i][k]*N[k][j];
    }
  }
  return L;
}

bool integer_coef(mi & A, int d){ //bool function returns true if 
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      if (A[i][j]%d == 0) A[i][j] = A[i][j]/d;
      else return false;
    }
  }
  return true;
}

void pushback(const mi & A){
  int i = 0;
  while (used[i]) i++;
  used[i] = 1;
  list[i] = A;
}

mi perm(const mi & N, long long n){
  mi M(4, vi(4));
  for (int i = 0; i < 4; i++){
    for(int j = 0; j < 4; j++){
      M[i][P[n][j]] = N[i][j];
    }
  }
  return M;
}

void transpose(mi & M){
  for (int i = 0; i < 4; i++){
    for(int j = 0; j < i; j++){
      long long a;
      a = M[i][j];
      M[i][j] = M[j][i];
      M[j][i] = a;
    }
  }
}

bool in_list(mi & A){
    for (int i = 0; i < 24 and used[i]; i++){
        if (A == list[i]) return true;
    }
    return false; 
}

int main(){
  initialize_P();
  string s;
  int iteration = 1;
  while (fe >> s){
    fs << "The matrices that send the first matrice into the second in case " << iteration << " are the following:" << endl << endl;
    ++iteration;
    list = vector<mi>(24, mi(4, vi(4,0)));
    used = vector<bool>(24,0);
    mi M(4, vi(4));
    read_matrix(M,s);
    long long dm = determinant(M);
    if (dm == 0) fs << "Your first matrix does not define a tetrahedron" << endl;  
    mi N(4, vi(4));
    fe >> s;
    read_matrix(N,s); 
    long long dn = determinant(M);
    if (dn == 0) fs << "Your second matrix does not define a tetrahedron" << endl;
  
    if (dn != dm and dn != -dm) fs << "They are not projectively equivalent" << endl;
  
    else{
      //we want A st A*M = N ===> A = N*M^(-1)
      mi Mai = almost_inverse_matrix(M);
      for(int i = 0; i < 24; i++){
        mi Naux = perm(N, i);
        mi A = multiply(Naux,Mai);
        if (integer_coef(A, dm)){
          if (not in_list(A)) pushback(A);
        }
      }
      int i = 0;
      while(used[i]){
        write_matrix(list[i]);
        fs << endl;
        i++;
      }
      if (i == 0) fs << "They are not projectively equivalent" << endl;
    }
    fe >> s;
    mi B(4, vi(4));
    read_matrix(B,s);
    transpose(B);
    bool b = 0;
  
    for(int i = 0; i < 24 and used[i]; i++){
      if (B == list[i]) b = 1;
    }
  
    if (b) cout << s << " is one of the transformations that send the first tetrahedron to the second" << endl;
    else cout << s <<  " is NOT one of the transformations that send the first tetrahedron to the second" << endl;
  }
  fe.close();
  fs.close();
}
