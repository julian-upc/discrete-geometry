#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <set>
#include <queue>

using namespace std;

typedef vector<double> V;//we must handle ring elements? only necessary for H4 
typedef vector<V> M;

/*
void escriuvector (V v){//writes a vector
	for (int i=0; i<v.size();++i){
		cout<<v[i]<<" ";
	}
	cout<<endl;
}
*/

struct veccomp{//defines how to compare two vectors in order to keep them ordered in the set
	bool operator() (const V& v1, const V& v2) const {
		int n=v1.size();
		for (int i=0; i<n; ++i){
			if (abs(v1[i]-v2[i])>0.00001){//hope
				return (v1[i]<v2[i]);
			}
		}
		return 0;
	}
};

M computematrixA (int &n){
	M An (n, V (n+2, 0));
	for (int i=0; i<n; ++i){
		An[i][i+1]=1;
		An[i][i+2]=-1;

	}
	return An;
}

M computematrixD (int &n){
	M Dn (n, V (n+1, 0));
	for (int i=0; i<n-1; ++i){
		Dn[i][i+1]=1;
		Dn[i][i+2]=-1;

	}
	Dn[n-1][n-1]=1;
	Dn[n-1][n]=1;
	return Dn;
}

M computematrixE (){
	M E (8, V (9, 0));
	for (int i=1; i<7; ++i){
		E[i-1][i]=1;
		E[i-1][i+1]=-1;
	}
	E[6][6]=1;
	E[6][7]=1;
	for (int i=1; i<9; ++i){
		E[7][i]=1;
	}
	return E;
}

M computematrixF (){
	M F (4, V (5, 0));
	F[0][1]=1;
	F[0][2]=-1;
	F[1][2]=1;
	F[1][3]=-1;
	F[2][3]=1;
	for (int i=1; i<5; ++i){
		F[3][i]=-1;
	}
	return F;
}

M computematrixH3 (){
	const double golden_ratio = 1.6180339887498948482;//taken from tables
	M H (3, V(4,0));
	H[0][1]=2;
	H[1][1]=-golden_ratio-1;
	H[1][2]=1;
	H[1][3]=-golden_ratio;
	H[2][3]=2;
	return H;
}

M computematrixH4 (){
	const double golden_ratio = 1.6180339887498948482;//taken from tables
	M H (4, V(5,0));
	H[0][1]=-golden_ratio-1;
	H[0][2]=1;
	H[0][3]=1;
	H[0][4]=1;
	for (int i=1; i<4; ++i){
		H[i][i]=-1;
		H[i][i+1]=1;
	}
	return H;
}


M computematrix (char& X, int& n){// computes the required matrix, if the entry is bad returns a zero value inside a matrix
	M zero (1,V(1,0));
	if (X=='A'){
		return computematrixA(n);
	}
	else if (X=='D'){
		return computematrixD(n);
	}
	else if (X=='E'){
		if (n==8) return computematrixE();
		else return zero;
	}
	else if (X=='F') {
		if (n==4) return computematrixF ();
		else return zero;
	}
	else if (X=='H'){
		if (n==3) return computematrixH3();
		else if (n==4) return computematrixH4();
		else return zero;
	}
	else return zero;
}

double scalarprod( V &v1,V& v2){//computes scalar product of two vectors
	int n= v1.size();
	double product=0;
	for (int i=0; i<n; ++i){
		product+=v1[i]*v2[i];
	}
	return product;
}

V multiply (double a, V v){//multiplies vector by scalar
	int n= v.size();
	for (int i=0; i<n; ++i){
		v[i]*=a;
	}
	return v; 
}

/*V choose (M &Xn){//chooses a random vector and checks if it lies in some hiperplane (the odds are bad so probably it won't run many times)
	int n=Xn[0].size();
	int m=Xn.size();
	V v (1,1);
	for (int i=1; i<n; ++i){
		double x = ((double) rand() / (RAND_MAX));
		double a= (double) i;
		v.push_back( x+a);//components of the vector random between i and i+1
	}
	for (int i=0; i<m; ++i){
		if (abs(scalarprod(v, Xn[i]))<0.1){//to avoid floating point trouble a tolerance margin would be good
			return choose (Xn);
		}
	}
	return v;
}*/

V choose (M &Xn){//aparently this vector is harmless for all matrices given except for F
	int n=Xn[0].size();
	V v (n,0);
	for (int i=0; i<n; ++i){
		v[i]= (double) i+1;
	}
	return v;
}

V chooseF (M& Xn){//chooses random vector for F
	V v(1,1);
	for (int i=1; i<5; ++i){
		double x = ((double) rand() / (RAND_MAX));
		double a= (double) i;
		v.push_back( x+a);//components of the vector random between i and i+1
	}
	for (int i=0; i<4; ++i){
		if (abs(scalarprod(v, Xn[i]))==0){//to avoid floating point trouble a tolerance margin would be good (not necessary for F apparently)
			v= chooseF (Xn);
		}
	}
	return v;
}

V substract( V v1, V v2){//subtracts two vectors (floating point trouble)
	int n=v1.size();
	V v; 
	for( int i=0; i<n; ++i){
		v.push_back(v1[i]-v2[i]);
	}
	return v;
}

V reflection (V vector, V hiperplane){//reflects the vector through hiperplane (begging for floating point trouble, look for simbolic calculations)
	int n=vector.size();
	V reflected_vector=vector ;
	double coeficient=2*scalarprod(vector,hiperplane)/scalarprod(hiperplane,hiperplane);
	hiperplane=multiply(coeficient, hiperplane);
	reflected_vector=substract(reflected_vector,hiperplane);
	return reflected_vector;
}

bool compare (V& vector, set<V, veccomp>& vectorlist){//is the vector in the vectorlist?
	set<V>::iterator it;
	it=vectorlist.find(vector);
	if (it==vectorlist.end()){
		return 0;
	}
	return 1;
} 

void  orbit (M &Xn, set<V, veccomp>& orbitlist, queue<int>& last_hyperplane_reflected, queue<V>& uncheckedorbit){//computes the orbit
	while (!uncheckedorbit.empty()){
		int n= Xn.size();
		for (int j=0; j<n; ++j){
			if (j!=last_hyperplane_reflected.front()){
				V newelement=reflection(uncheckedorbit.front(),Xn[j]);
				if (compare(newelement,orbitlist)==0){
					orbitlist.insert(newelement);
					uncheckedorbit.push(newelement);
					last_hyperplane_reflected.push(j);
				}
			}
		}
		uncheckedorbit.pop();
		last_hyperplane_reflected.pop();
	}
}

/*
void escriu (M X){ // writes a matrix
	cout<<"La matriu es "<<endl;
	for (int i=0; i<X.size(); ++i){
		for (int j=0; j<X[0].size(); ++j){
			cout<<" "<<X[i][j];
		}
		cout<<endl;
	}

}
*/

int main(){
	srand (time(NULL));//random seed in case random numbers are necessary such as chooseF
	char X;
	int n;
	double firstcoordinate;
	V v;
	cin>>X>>n>>firstcoordinate;
	M zero (1,V(1,0));
	M Xn=computematrix(X,n);
	
	if (Xn.size()!=1){
		//escriu(Xn);
		if(firstcoordinate!=0){
			v.push_back(firstcoordinate);
			for(int i=0; i<(Xn[0].size()-1); ++i){
				double a;
				cin>>a;
				v.push_back(a);
			}
		}
		else{
			if (X=='F') v=chooseF(Xn);
			else v=choose(Xn);
			//escriuvector (v);
		}
		set<V, veccomp> orbitlist;
		queue <V> uncheckedorbit;
		orbitlist.insert(v);
		uncheckedorbit.push(v);
		queue <int> last_hyperplane_reflected;
		last_hyperplane_reflected.push(-1);
		orbit (Xn, orbitlist, last_hyperplane_reflected, uncheckedorbit);
		cout<<"Orbit size is "<<orbitlist.size()<<endl;
	}
	else{
		cout<<"This is not a valid input"<<endl; // error message 
	}
	
}
