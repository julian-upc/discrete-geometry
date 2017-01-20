//
//  intersectionsTest.cpp
//  PolytopeProject
//


#include "intersectionsTest.hpp"

Intersections::Intersections():p(Point(V{0,0}, true)),q(Point(V{2,0}, false)){
	
	M = vector<VV>(6, VV(6, V{0,0}));
	n_gale_diagrams = VVI(20, VI(2,0));

	
	VP v(6, Point(V{0,0}, true));
	
	for(int i = 0; i<6; ++i) v[i] = Point(V{1,i+1}, true);
	
	VP l(6, Point(V{0,0}, true));
	VP r(6, Point(V{0,0}, true));
	
	V p1 = p.get_coordinates();
	V q1 = q.get_coordinates();
	
	for(int i = 0; i < 6; ++i){
		V p2 = v[i].get_coordinates();
		V q2 = v[i].get_coordinates();
		
		mpq_class m1 = (p2[1]-p1[1])/(p2[0]-p1[0]);
		mpq_class m2 = (q2[1]-q1[1])/(q2[0]-q1[0]);
		
		m1.canonicalize();
		m2.canonicalize();
		
		mpq_class c1 = m1*p1[0];
		mpq_class c2 = m2*q1[0];
		
		c1.canonicalize();
		c2.canonicalize();
		
		c1 = p1[1]-c1;
		c2 = q1[1]-c2;
		
		c1.canonicalize();
		c2.canonicalize();
		
		l[i] = Point(V{m1,c1}, true);
		r[i] = Point(V{m2,c2}, true);
	}
	
	for(int i = 0; i < 6; ++i){
		for(int j = 0; j < 6; ++j){	
			V r1 = l[i].get_coordinates();
			V r2 = r[j].get_coordinates();
			
			mpq_class x1 = -r2[1]+r1[1];
			mpq_class x2 = r2[0]-r1[0];
			
			x1.canonicalize();
			x2.canonicalize();
			
			mpq_class x = x1/x2;
			x.canonicalize();
			
			mpq_class y = r1[0]*x;
			y.canonicalize();
			
			y = y + r1[1];
			y.canonicalize();
			
			M[i][j] = V{x,y};
		}
	}
	
	VI vi(6,0);
	VB b(6,true);
	VB vs(8, false);
	
	p_combinations(vi, b, p, q, 6);
	sign_combinations(vs, 7, 5, 4);
}

void Intersections::sign_combinations(VB& v, int n, int ps, int ns){
    if(n>0){
        if(ps>0){
            v[n-1] = true;
            sign_combinations(v, n-1, ps-1, ns);
        }
        if(ns>0){
            v[n-1] = false;
            sign_combinations(v, n-1, ps, ns-1);
        }
    }
    else{
        all_s_combinations.push_back(VB(v));
    }
}

void Intersections::p_combinations(VI& v, VB& b, const Point& p, const Point& q, int n){
    if(n>0){
        for(int i = 0; i < 6; ++i){
            if(b[i]){
                v[n-1] = i;
                b[i]=false;
                p_combinations(v, b, p, q, n-1);
                b[i]=true;
            }
        }
    }
    else all_p_combinations.push_back(VI(v));
}

VGD Intersections::neighborly_diagrams_hardcoded(){
	
	n_gale_diagrams[0] = {0, 34};
	n_gale_diagrams[1] = {5, 34};
	n_gale_diagrams[2] = {16, 34};
	n_gale_diagrams[3] = {103, 34};
	n_gale_diagrams[4] = {119, 34};
	n_gale_diagrams[5] = {138, 84};
	n_gale_diagrams[6] = {143, 84};
	n_gale_diagrams[7] = {264, 34};
	n_gale_diagrams[8] = {269, 34};
	n_gale_diagrams[9] = {288, 34};
	n_gale_diagrams[10] = {304, 34};
	n_gale_diagrams[11] = {337, 34};
	n_gale_diagrams[12] = {361, 84};
	n_gale_diagrams[13] = {380, 84};
	n_gale_diagrams[14] = {415, 84};
	n_gale_diagrams[15] = {506, 34};
	n_gale_diagrams[16] = {508, 34};
	n_gale_diagrams[17] = {566, 34};
	n_gale_diagrams[18] = {576, 34};
	n_gale_diagrams[19] = {600, 84};
	
	VGD vgd;
	
	for(int i = 0; i < n_gale_diagrams.size(); ++i){
		
		int n = n_gale_diagrams[i][0];
		int m = n_gale_diagrams[i][1];
		
		VI pcomb = all_p_combinations[n];
		VB scomb = all_s_combinations[m];
		
		VP vp(8, Point(V{0,0}, false));
		
		for (int k=0; k < pcomb.size(); ++k) vp[k] = Point(M[k][pcomb[k]], scomb[k]);
		
		vp[6] = Point(V{0,0}, scomb[6]);
		vp[7] = Point(V{2,0}, scomb[7]);
		
		Gale_Diagram gd(vp);
		
		gd.print();
		
		vgd.push_back(gd);
	}
	
	return vgd;
}

VGD Intersections::neighborly_diagrams(){
	
	VP vp(8, Point(V{0,0}, false));
	vp[6] = p;
	vp[7] = q;
	
	VGD vgd;
	
	for(int i = 0; i < all_p_combinations.size(); ++i){
		for(int j = 0; j < all_s_combinations.size(); ++j){
				VI p_combination = all_p_combinations[i];
				VB s_combination = all_s_combinations[j];
				for(int k = 0; k < 6; ++k) vp[k] = Point(M[k][p_combination[k]],s_combination[k]);
				vp[6].set_sign(s_combination[6]);
				vp[7].set_sign(s_combination[7]);
				
				++number_diagrams;
				
				Gale_Diagram gd(vp);
				
				if(gd.isSimplicial()){
					++number_simplicial;
					if(gd.is_neighborly()){
						++number_neighborly;
						vgd.push_back(gd);
					}
				}
		}
	}
	
	cout << "totals: " << endl;
	cout << number_diagrams << " " << number_simplicial << " " << number_neighborly << endl;
	
	return vgd;
}




