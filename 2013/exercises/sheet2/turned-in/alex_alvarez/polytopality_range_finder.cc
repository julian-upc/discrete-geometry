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
typedef adjacency_list<vecS, vecS, undirectedS, property<vertex_index_t, int> > boost_graph;

//Adjacency list of the graph
VVI adj;

//Checks planarity of the graph using Boost's implementation of Boyer-Myrvold's algorithm
bool check_planarity() {
  int n = adj.size();
  boost_graph G(n);
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
  for (int i = 0; i + 1 < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      int dif = j - i, dif2 = i - j + n;
      bool valid_edge = false;
      for (int k = 0; k < int(S.size()) and not valid_edge; ++k)
        if (S[k] == dif or S[k] == dif2) valid_edge = true;
      if (valid_edge) {
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

//Return the number of bits in x set to one
int set_bits(int x) {
  int res = 0;
  while (x > 0) {
    if (x & 1) ++res;
    x >>= 1;
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
    if (set_bits(mask) == d - 1) {
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
//Backtracks to find a principal subdivision of the complete graph
bool rec(int act, int u, set<int>& used, vector<PII>& connections) {
  if (u == connections[act].second) {
    if (act + 1 == connections.size()) return true;
    return rec(act + 1, connections[act + 1].first, used, connections);
  }
  bool found = false;
  for (int i = 0; i < adj[u].size() and not found; ++i) {
    int v = adj[u][i];
    if (used.count(v) and connections[act].second != v) continue;
    used.insert(v);

    if (rec(act, v, used, connections)) found = true;
    if (connections[act].second != v) used.erase(v);
  }
  return found;
}

//Checks with a backtracking if there is a principal subdivision of K_{d + 1}
bool psp(int d) {
  if (adj[0].size() < d) return false;
  vector<int> vertices;
  for (int i = 0; i < d; ++i) vertices.push_back(adj[0][i]);
  set<int> used;
  
  used.insert(0);
  used.insert(adj[0].begin(), adj[0].begin() + d);
  vector<PII> connections;
  for (int i = 0; i + 1 < d; ++i) {
    for (int j = i + 1; j < d; ++j) {
      bool already_connected = false;
      for (int k = 0; k < adj[vertices[i]].size(); ++k) {
        if (adj[vertices[i]][k] == vertices[j]) {
          already_connected = true;          
        }
      }
      if (not already_connected) connections.push_back(PII(vertices[i], vertices[j]));
    }
  }
  
  return connections.size() == 0 ? true : rec(0, connections[0].first, used, connections);
}


int main(int argc, char** argv) {
  //Number of nodes
  int n;
  
  //Size of S
  int s_size;
  
  if (argc < 3 or argc < atoi(argv[2]) + 3) {
    cout << "usage: program n k s1 ... sk" << endl;
    cout << "\t\tn: the number of nodes" << endl;
    cout << "\t\tk: the number of differences that define the circulant graph" << endl;
    cout << "\t\tsi: the i-th difference" << endl;
    return 1;
  }
  
  n = atoi(argv[1]);
  s_size = atoi(argv[2]);
  
  VI S(s_size);
  for (int i = 0; i < s_size; ++i) 
    S[i] = atoi(argv[2 + i]);
  
  //Generate the graph
  generate_graph(n, S);
  
  for (int i = 0; i < n; ++i) {
    cerr << i << ": ";
    for (int j = 0; j < adj[i].size(); ++j) cerr << adj[i][j] << " ";
    cerr << endl;
  }
  
  //candidate[i] = true iff i is a valid dimension so far
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