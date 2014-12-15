#include <iostream>
#include <set>
#include <queue>
#include "reflection_group.cc"

using namespace std;

Vec reflect(const Vec& v, const Vec& h) {
  Base h2 = 0;
  for (int i=0; i<h.size(); i++) {
    h2 += h[i]*h[i];
  }
  
  Base hv = 0;
  for (int i=0; i<v.size(); i++) {
    hv += h[i]*v[i];
  }
  
  Vec reflexion(v.size());
  for (int i=0; i<v.size(); i++) {
    reflexion[i] = v[i] - 2*hv*h[i]/h2;
  }
  
  return reflexion;
}

set<Vec> orbit(const Mat& m, const Vec& v) {
  set<Vec> s;
  queue<Vec> q;
  s.insert(v);
  q.push(v);

  while (!q.empty()) {
    Vec current_v = q.front();
    q.pop();
    
    for(int i=0; i<m.size(); i++) {
      Vec reflexion = reflect(current_v, m[i]);
      set<Vec>::iterator it = s.find(reflexion);
      if (it == s.end()) {
        s.insert(reflexion);
        q.push(reflexion);
    
        if (s.size() % 100000 == 0) {
          cout << "Current set size: " << s.size() << endl;
        }
      }
    }
  }

  return s;
}

int main() {
  char c;
  int n;
  
  cin >> c >> n;
  
  if (n >= 1 && (c == 'A' || c == 'D' || (c == 'E' && n == 8) || (c == 'F' && n == 4) || (c == 'H' && (n == 3 || n == 4)))) {
    Mat m = generators(c, n);
    show(m);
    Vec v(m[0].size());

    if (cin >> v[0]) {
      for (int i=1; i<v.size(); i++) {
        cin >> v[i];
      }
    } else {
      v[0] = 1;
      for (int i=1; i<v.size(); i++) {
        v[i] = i;
        //v[i] = 0;
        //for (int j=0; j<m.size(); j++) {
        //  v[i] += m[j][i];
        //}
      }
    }
    
    set<Vec> s = orbit(m, v);
    cout << "|" << c << n << "| = " << s.size() << endl;
  } else {
    cout << "Incorrect input" << endl;
    return 0;
  } 
}
