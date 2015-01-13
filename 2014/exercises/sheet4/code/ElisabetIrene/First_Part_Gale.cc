#include <iostream>
#include <math.h>
#include <vector> 
#include <set> 
#include <algorithm> 
#include <stdexcept>

//using namespace std;

typedef std::vector< std::vector< std::vector<int> > > Tri_vector;
typedef std::vector< std::vector<int> >  Matriz;

bool are_different( std::vector<int> v1, std::vector<int> v2) { 
    int minim;
    minim = std::min(v1.size(), v2.size());
    for (int i = 0; i < minim; ++i) {
        if (v1[i] != v2[i]) return true;
    }
    return false;
}

int compare( int a, int b, int c ) {
    if ( a <= 0 || b <= 0 || c < 0 ) {
        throw std::invalid_argument( "Received invalid values" );
    }
}

void write_tri(const Tri_vector& t) {
  for (int i = 0; i < t.size(); ++i) {
    for (int j = 0; j < t[i].size(); ++j) {
      for (int k = 0; k < t[i][j].size(); ++k) {
        if (k != 0) std::cout << ' ';
        std::cout << t[i][j][k];
      }
      std::cout << std::endl;
    }
    std::cout << "-----" << std::endl;
  }
}

void write_mat(const Matriz& m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m[i].size(); ++j) {
            if (j != 0) std::cout << ' ';
            std::cout << m[i][j];
        }
        std::cout << std::endl;
    }
} 

void write_mat(const std::vector< std::vector<bool> >& m) {
    for (int i = 0; i < m.size(); ++i) {
        for (int j = 0; j < m[i].size(); ++j) {
            if (j != 0) std::cout << ' ';
            std::cout << m[i][j];
        }
        std::cout << std::endl;
    }
}

void write_vec(const std::vector<int>& v) {
  for (int i = 0; i < v.size(); ++i) {
    if (i != 0) std::cout << ' ';
    std::cout << v[i];
  }
  std::cout << std::endl;
}

int main() {
  int e, n, m, kmax, gcompl, val_actual, val_0, Gsize;
  std::set< std::vector<int> > Coor1D;
  std::vector<int> coor_actual;
  std::vector<bool> lexi_actual;
  std::vector< std::vector<int> > conf_actual, conf_0;
  std::vector< std::vector<bool> > Gale_lexi, Gale_lexi0;
  Tri_vector Sk_1, Sk, Gale_conf, Gale_0;
  bool is_lexi;
  // comprobar entrada: enteros positivos
  std::cout << "Introduce e, n y m:" << std::endl;
  std::cin >> e >> n >> m;
  compare(e, n, m);
  kmax = floor(n/2)*m;
  // Generación de configuraciones unidimensionales
  // inicializacion
  if (n > 1) {
    Sk_1.resize(n-1);
    Sk.resize(n-1);
    for (int g = 0; g < n-1; ++g) {
      Sk_1[g].resize(1);
      Sk_1[g][0].resize(g+1);
    }
    // emparejado primero
    for (int g = 0; g < 1; ++g) {
      gcompl = n - 2 - g;
      coor_actual = Sk_1[g][0];
      coor_actual.insert(coor_actual.end(),Sk_1[gcompl][0].begin(),Sk_1[gcompl][0].end()); 
      Coor1D.insert(coor_actual);
    }
  }
  else {
    coor_actual.clear();
    coor_actual.resize(1);
    Coor1D.insert(coor_actual);
  }
  //write_mat(Coor1D);



  for (int k = 1; k < kmax + 1; ++k) {
    for (int g = 0; g < n-1; ++g) { // el nuemro de puntos en la config
      for (int i = 0; i < Sk_1[g].size(); ++i) { // el numero de opciones 
	val_0 = Sk_1[g][i][0];	
	for (int j = 0; j < g + 1; ++j) { // las coordenadas, la que aumentar
	  val_actual = Sk_1[g][i][j];
	  if ( (j == 0) && (val_actual < m) ) {
	    coor_actual = Sk_1[g][i];
	    coor_actual[0] += 1;
	    Sk[g].push_back(coor_actual);
          }
	  else if ( ((val_0 - val_actual) == 1) && (val_actual < m) && ((Sk_1[g][i][j-1] - val_actual) == 1) ){
	    coor_actual = Sk_1[g][i];
	    coor_actual[j] += 1;
	    Sk[g].push_back(coor_actual);
          }
	}
      }
    }
   // emparejado 
    for (int g = 0; g < ceil((n-1)/2); ++g) {
      gcompl = n - 2 - g;
      for (int i = 0; i < Sk[g].size(); ++i) {
        for (int j = 0; j < Sk[gcompl].size(); ++j) {
  	  coor_actual = Sk[g][i];
	  for (int l = 0; l < coor_actual.size(); ++l) {
	    coor_actual[l] *= -1;
	  }
	  coor_actual.insert(coor_actual.end(),Sk[gcompl][j].begin(),Sk[gcompl][j].end());
	  std::sort(coor_actual.rbegin(), coor_actual.rend()); // reverse iterators
	  //Coor1D.push_back(coor_actual);
	  Coor1D.insert(coor_actual);	  
	  coor_actual = Sk[gcompl][j];
	  if (are_different(coor_actual, Sk[g][i])) {
	    for (int l = 0; l < coor_actual.size(); ++l) {
	      coor_actual[l] *= -1;
	    }
	    coor_actual.insert(coor_actual.end(),Sk[g][i].begin(),Sk[g][i].end());
	    std::sort(coor_actual.rbegin(), coor_actual.rend()); // reverse iterators
	    Coor1D.insert(coor_actual);
	  }
        }
      }
    }
    Sk_1 = Sk; 
    Sk.clear();
    Sk.resize(n-1);
  }
  
  // Primera dimension
  conf_actual.resize(n);
  lexi_actual.resize(n);
  for (int i = 0; i < n; ++i) {
    conf_actual[i].resize(e);
  }
  for (std::set< std::vector<int> >::iterator it = Coor1D.begin(); it != Coor1D.end(); ++it) {
    coor_actual = (*it);
    if (coor_actual[0] == m) {
      for (int i = 0; i < n; ++i) {
   	conf_actual[i][0] = coor_actual[i];
	if (i == 0) {
	  lexi_actual[0] = true; 
	}      
	else if ((coor_actual[i-1] - coor_actual[i]) > 0) {
	  lexi_actual[i] = true;
	}
	else {
	  lexi_actual[i] = false;
	}
      }
      Gale_conf.push_back(conf_actual);
      Gale_lexi.push_back(lexi_actual);
    }   
    //conf_actual.clear();
  }
  // Ampliación hasta e dimensiones
  for (int d = 1; d < e; ++d) {
    Gale_0 = Gale_conf;
    Gale_lexi0 = Gale_lexi;
    Gale_conf.clear();
    Gale_lexi.clear();
    for (std::set< std::vector<int> >::iterator it = Coor1D.begin(); it != Coor1D.end(); ++it) {
      coor_actual = (*it);
      for (int i = 0; i < Gale_0.size(); ++i) {
        conf_0 = Gale_0[i];
        lexi_actual = Gale_lexi0[i];
        do {
 	  is_lexi = true;
          if (conf_0[0][d-1] < coor_actual[0]) {
	    is_lexi = false;
	    continue;
	  }
	  for (int j = 1; j < n; ++j) {
            if ( (!lexi_actual[j]) && ((coor_actual[j-1] - coor_actual[j]) < 0) ) {
	      is_lexi = false;
	      break;
	    }
	    else if ( (!lexi_actual[j]) && ((coor_actual[j-1] - coor_actual[j]) > 0) ) {
	      lexi_actual[j] = true;
  	    }
          }
	  if (is_lexi) {
	    conf_actual = conf_0;
	    for (int j = 0; j < n; ++j) {
	      conf_actual[j][d] = coor_actual[j];
	    }
	    //añadir conf_actual y lexi_actual
	    Gale_conf.push_back(conf_actual);
	    Gale_lexi.push_back(lexi_actual);
	  }
        } while(std::prev_permutation(coor_actual.begin(), coor_actual.end()));
      }
    }
  }
  write_tri(Gale_conf);

}
