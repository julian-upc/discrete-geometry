#include<iostream>
#include<vector>
#include"auxiliar.cpp"

using namespace std;

int reflexion_a (int a, int n,int i){ //reflexion with axis going through two vertices. i ranges from 0 to n/2-1
	int x=(a+i)%n;
	x=n-x;
	x=x-i;
	if (x<0) x+=n;
	return x;
}

int reflexion_b (int a, int n,int i){
	int x=(a+i)%n;
	x=n-1-x;
	x=x-i;
	if (x<0) x+=n;
	return x;
}


bool check_dihedral(vector<bool> v1, vector<bool> v2) {  //returns true if v1 and v2 are the same up to permutations in the diedral group
	int n=v1.size();
	for (int i=1; i<n; ++i){  // check rotations of angle 2*pi*i/n
		bool equivalent=true;
		for (int j=0; j<n; ++j){
			if (v1[j]!=v2[(j+i)%n]) {    //check if element j of the v1 and its image on v2 are the same
				equivalent=false;
				break;
			}
		}
		if (equivalent) return true;
	}
	
	for (int i=1; i<n/2; ++i){ // check reflexions
		bool equivalent=true;
		for (int j=0; j<n; ++j){  //check reflexions of type a
			if (v1[j]!=v2[reflexion_a(j,n,i)]){
				equivalent=false;
				break;
			}
		}
		if (equivalent) return true;
		equivalent=true;
		for(int j=0; j<n; ++j){
			if (v1[j]!=v2[reflexion_b(j,n,i)]){
				equivalent=false;
				break;
			}
		}
		if (equivalent) return true;
	}
	return false;
}

bool check_dihedral(vector<PointConfig> v1, vector<bool> u2) {  //returns true if v1 and v2 are the same up to permutations in the diedral group
	int n=v1.size();
	vector<bool> u1 (n);
  //, u2(n);
	for (int i=0; i<n; ++i){
		u1[i]=v1[i].positive;
		//u2[i]=v2[i].positive;
	}
	return check_dihedral(u1,u2);
}


/*int main(){
	vector<bool> v1{1,0,1,0,0,1,0,0};
	vector<bool> v2{0,0,1,0,0,1,0,1};
	if (check_dihedral(v1,v2)) cout << "same";
	else cout << "dif";
}*/


