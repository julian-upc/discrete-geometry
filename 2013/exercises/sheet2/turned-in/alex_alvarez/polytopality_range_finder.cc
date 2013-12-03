#include <iostream>
#include <vector>
#include <set>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>
using namespace std;
using namespace boost;

typedef vector<int> VI;
typedef vector<VI> VVI;
typedef pair<int, int> PII;
typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int> > graph_boost;

//Adjacency list of the graph
VVI adj;

//Checks planarity of the graph using Boost's implementation of Boyer-Myrvold's algorithm
bool check_planarity() {
  int n = adj.size();
  graph_boost G(n);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < int(adj[i].size()); ++j) {
      int v = adj[i][j];
      if (i < v) add_edge(i, v, G);
    }
  }
  return boyer_myrvold_planarity_test(G);
}

//Generates the graph of n vertices and S as subset of differences
void generate_graph(int n, const VI& S) {
  adj = VVI(n);
  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int dif = j - i, dif2 = i - j + n;
      bool ok = false;
      for (int k = 0; k < int(S.size()) and not ok; ++k)
        if (S[k] == dif or S[k] == dif2) ok = true;
      if (ok) {
        adj[i].push_back(j);
        adj[j].push_back(i);
      }
    }
  }
}

//A simple DFS to count the reachable nodes
int dfs(int u, vector<bool>& seen, int mask) {
  seen[u] = true;
  int res = 1;
  for (int i = 0; i < int(adj[u].size()); ++i) {
    int v = adj[u][i];
    if (not seen[v] and ((mask >> v)&1) == 0) res += dfs(v, seen, mask);
  }
  return res;
}

//Checks Balinski's theorem for dimension d
bool check_balinski(int d) {
  int n = adj.size();
  
  //Iterate (in a very inefficient way but fast to code way)
  //over every (d-1)-subset of the nodes
  for (int mask = 0; mask < (1 << n); ++mask) {
    //The bits sets to one codify the cut-set
    if (__builtin_popcount(mask) == d - 1) {
      vector<bool> seen(n, false);
      
      //Select a vertex that is not in the candidate cut
      int start_vertex = 0;
      while ((mask >> start_vertex)&1) ++start_vertex;
      
      //DFS to check whether there's more than one component
      if (dfs(start_vertex, seen, mask) != n - d + 1) return false;
    }
  }
  return true;
}

//Auxiliar recursive function for the PSP checker
bool rec(int act, int u, set<int>& usat, vector<PII>& connections) {
  if (u == connections[act].second) {
    if (act + 1 == connections.size()) return true;
    return rec(act + 1, connections[act + 1].first, usat, connections);
  }
  bool ok = false;
  for (int i = 0; i < adj[u].size() and not ok; ++i) {
    int v = adj[u][i];
    if (usat.count(v) and connections[act].second != v) continue;
    usat.insert(v);

    if (rec(act, v, usat, connections)) ok = true;
    if (connections[act].second != v) usat.erase(v);
  }
  return ok;
}

//Checks with a backtracking if there is a principal subdivision of K_{d + 1}
bool psp(int d) {
  if (adj[0].size() < d) return false;
  vector<int> vertices;
  for (int i = 0; i < d; ++i) vertices.push_back(adj[0][i]);
  set<int> usat;
  
  usat.insert(0);
  for (int i = 0; i < d; ++i) usat.insert(adj[0][i]);
  vector<PII> connections;
  for (int i = 0; i < d; ++i) {
    for (int j = i + 1; j < d; ++j) {
      bool ok = false;
      for (int k = 0; k < adj[vertices[i]].size(); ++k) {
        if (adj[vertices[i]][k] == vertices[j]) {
          ok = true;          
        }
      }
      if (not ok) connections.push_back(PII(vertices[i], vertices[j]));
    }
  }
  
  if (connections.size() == 0) return true;
  return rec(0, connections[0].first, usat, connections);
}


int main() {
  //Number of nodes
  int n;
  cin >> n;
  
  //Size of S
  int s_size;
  cin >> s_size;
  
  VI S(s_size);
  for (int i = 0; i < s_size; ++i) 
    cin >> S[i];
  
  //Generate the graph
  generate_graph(n, S);
  
  for (int i = 0; i < n; ++i) {
    cerr << i << ": ";
    for (int j = 0; j < adj[i].size(); ++j) cerr << adj[i][j] << " ";
    cerr << endl;
  }
  
  vector<bool> candidate(n, true);
  candidate[0] = candidate[1] = false;
  
  //Check the number of edges for d = 2
  int edges = 0;
  for (int i = 0; i < n; ++i) edges += adj[i].size();
  edges /= 2;
  if (n - edges != 0) candidate[2] = false;  
  
  //Check in a naive way Balinski theorem
  for (int i = 2; i < n; ++i) {
    if (candidate[i]) candidate[i] = check_balinski(i);    
  }
  
  for (int i = 4; i < n; ++i) if (candidate[i]) candidate[i] = psp(i);
  
  if (candidate[3]) candidate[3] = check_planarity();
    
  cout << "Valid candidates so far:";
  for (int i = 0; i < n; ++i) {
    if (candidate[i]) cout << " " << i;
  }
  cout << endl;
}