//
//  Gale_Diagram.cpp
//  PolytopeProject
//


#include "Gale_Diagram.hpp"



Gale_Diagram::Gale_Diagram(vector<Point> points): points(points),matrix(Matrix(points)){
};

bool Gale_Diagram::is_neighborly() const{
    //Choose 2 points
    for(int i=0; i<points.size()-1; i++){
        for(int j=i+1; j<points.size(); j++){
            
            bool edgeIJPresent=false;
            
            
            //Then choose another 2
            for(int k=0; k<points.size()-1 && (!edgeIJPresent); k++){
                for(int l=k+1; l<points.size() && (!edgeIJPresent); l++){
					if( (i!=k) and (i!=l) and (j!=k) and (j!=l)){
						edgeIJPresent = isSimplex(i, j, k, l);
					}
                }
            }
            
            if(!edgeIJPresent){
                return false;
            }
        }
        
    }
    return true;
}



vector<int> Gale_Diagram::findMissing(int i, int j, int k, int l) const{
    vector<int> pointNums(points.size());
    for(int a=0; a<points.size(); a++) pointNums[a]=a;

	pointNums[i] = -1;
	pointNums[j] = -1;
	pointNums[k] = -1;
	pointNums[l] = -1;
	
    vector<int> pN;
    
    for(int a=0; a<points.size(); a++){
		if(pointNums[a]>=0) pN.push_back(a);
	}
    
    return pN;
}


//Does NOT check for degeneracies
//returns -1 if x lies on the right of the oriented line ab, otherwise +1;
int Gale_Diagram::orientation(const Point &x, const Point &a, const Point &b) const{
    Matrix M = Matrix(vector<Point>{b.minus(a),x.minus(a)});
    return sgn(M.determinant());
}
//Checks if onePoint lies in the triangle defined by threePoints
bool Gale_Diagram::in_triangle(const vector<Point> &onePoint, const vector<Point> &triangle) const{
    Point x=onePoint[0];
    Point a=triangle[0];
    Point b=triangle[1];
    Point c=triangle[2];
    int o1 = orientation(x,a,b);
    int o2 = orientation(x,b,c);
    int o3 = orientation(x, c,a);
    return o1==o2 && o2 == o3;
}
// returns whether both linesegments (vectors of size 2) intersect.
bool Gale_Diagram::intersectingLineSegments(const vector<Point> &lineSegment1, const vector<Point> &lineSegment2) const{
    Point a=lineSegment1[0];
    Point b=lineSegment1[1];
    Point c=lineSegment2[0];
    Point d=lineSegment2[1];
    int o1 = orientation(a, c, d);
    int o2 = orientation(b,c,d);
    int o3 = orientation(c,a,b);
    int o4 = orientation(d,a,b);
    return o1!=o2 && o3!=o4;
}
bool Gale_Diagram::check_intersecting(const vector<Point> &points) const{
    if(points.size()!=4) throw runtime_error("Wrong Size in check_intersecting, Gale_diagram");
    
    vector<Point> positives;
    vector<Point> negatives;
    for(Point p : points){
        if(p.get_sign()) positives.push_back(p);
        else negatives.push_back(p);
    }
    int s = (int) positives.size();
    if(s==0 || s==4) return false;
    if(s==1) return in_triangle(positives, negatives);
    if(s==2) return intersectingLineSegments(positives, negatives);
    return in_triangle(negatives, positives);
    return true;
}

void Gale_Diagram::test() const{
    Point a(vector<double>{0,0},true);
    Point b(vector<double>{0,1},true);
    Point c(vector<double>{1,0},true);
    Point d(vector<double>{0,-1},true);
    Point e(vector<double>{-1,0},true);
    Point f(vector<double>{1,1},true);
    
    cout << (orientation(b,a,c)==1) << (orientation(d,a,c)==-1) << (orientation(c,a,b)==-1);
    cout << (orientation(e,a,b)==1) << endl;
    
    vector<Point> zero{a};
    vector<Point> triangle{b,f,c};
    vector<Point> triangle2{f,d,e};

    cout << "triangle results " << (in_triangle(zero, triangle)==false) << (in_triangle(zero, triangle2)==true) << endl;
    
}


bool Gale_Diagram::isSimplex(int i, int j , int k, int l) const{
    //Remove the columns i,j,k,l by selecting the other columns:
    vector<int> pointIndices = findMissing(i,j,k,l);
    /*cout << i << " " << j << " " << k << " " << l << endl;
    for(int m=0; m<pointIndices.size(); m++){
		cout << pointIndices[m] << " ";
	}
	cout << endl;*/
    //Select points
    vector<Point> fourPoints;
    for(int m=0; m<pointIndices.size(); m++){
        int pointIndex = pointIndices[m];
        Point p(points[pointIndex].get_coordinates(), points[pointIndex].get_sign());
        fourPoints.push_back(p);
    }
    
    return check_intersecting(fourPoints);
}
void Gale_Diagram::writeGraph(string filename){
    FILE * pFile;
    
    pFile = fopen((filename+".txt").c_str() , "w");
    igraph_t graph = makeVertexFacetStructure();
    igraph_write_graph_edgelist(&graph,pFile);
    igraph_destroy(&graph);
    fclose (pFile);
}

igraph_t Gale_Diagram::makeVertexFacetStructure() const{
    int s = (int)points.size();
    
    igraph_t graph;
    igraph_empty(&graph, 28, IGRAPH_UNDIRECTED);
    
    int facetCounter=s;
    int counter = 0;
    
    for(int i=0; i<s-3; i++){
        for(int j=i+1; j<s-2; j++){
            for(int k=j+1; k<s-1; k++){
                for(int l=k+1; l<s; l++){
                    counter++;

                    if(isSimplex(i, j, k, l)){
                        igraph_vector_t edges;
                        
                        igraph_vector_init(&edges, 8);
                        
                        vector<int> pointIndicesSimplex = findMissing(i,j,k,l);
                        for(int i = 0; i<4; i++){
                            VECTOR(edges)[2*i] = facetCounter;
                            VECTOR(edges)[2*i+1] = pointIndicesSimplex[i];
                        }
                        facetCounter++;
                        
                        igraph_add_edges(&graph, &edges, 0);
                        
                        //writeGraph(graph, ("counting"+to_string(counter)+".txt").c_str());
                        
                        igraph_vector_destroy(&edges);
                    }
                }
            }
        }
    }
    
    igraph_to_undirected(&graph, IGRAPH_TO_UNDIRECTED_COLLAPSE, 0);
    
    
    return graph;
}

bool Gale_Diagram::isIsomorphic(const Gale_Diagram &gale) const{
    igraph_bool_t result = false;
    igraph_t graph = makeVertexFacetStructure();
    igraph_t otherGraph = gale.makeVertexFacetStructure();
    igraph_isomorphic(&graph, &otherGraph, &result);
    return result;
}

//For any three points, looks whether they are on a line.
bool Gale_Diagram::isSimplicial() const{
    int s = (int)points.size();
    for(int i=0; i<s-3; i++){
        for(int j=i+1; j<s-2; j++){
            //If ij is a vertical line
            if(points[i].get_coordinates()[0]==points[j].get_coordinates()[0]){
                for(int k=j+1; k<s-1; k++){
                    if(points[k].get_coordinates()[0]==points[i].get_coordinates()[0]) return false;
                }
                continue;
            }
            
            //check if k is on line ij
            Point slopeCoordinates = points[i].minus(points[j]);
            mpq_class slope = slopeCoordinates.get_coordinates()[1]/slopeCoordinates.get_coordinates()[0];
            slope.canonicalize();
            
            mpq_class aSlope = slope*points[i].get_coordinates()[0];
            aSlope.canonicalize();
            aSlope -= points[i].get_coordinates()[1];
            aSlope.canonicalize();
            
            for(int k=j+1; k<s-1; k++){
                mpq_class rightHandSide = slope*points[k].get_coordinates()[0];
                rightHandSide.canonicalize();
                
                rightHandSide-=aSlope;
                rightHandSide.canonicalize();
                
                if(points[k].get_coordinates()[1] == rightHandSide){
                    return false;
                }
            }
        }
    }
    return true;
}




Matrix Gale_Diagram::galeToPolytope() const{
    //Make points homogene: positive points on z=1 and negative points on z=-1
    vector<vector<mpq_class>> B;
    for (int i = 0; i < 8; ++i){
        vector<mpq_class> homogeneCoordinates = points[i].get_coordinates();
        if(points[i].get_sign()) {
            homogeneCoordinates.push_back(1);
        }
        else{
            homogeneCoordinates.push_back(-1);
        }
        B.push_back(homogeneCoordinates);
    }
    
    Matrix BT = Matrix(B).transpose();
    vector<vector<mpq_class>> Btransposed=BT.getMatrix();
    
    
    //balance gale vector configuration
    vector<vector<mpq_class>> MPoints;
    for(int i = 0; i<3; ++i) {
        MPoints.push_back(B[i]);
    }
    
    Point fourthColumn(B[3],true);
    fourthColumn = fourthColumn.multiply(-1);
    for(int i = 4; i < 8; ++i){
        fourthColumn = fourthColumn.minus(Point(B[i],true));
    }
    MPoints.push_back(fourthColumn.get_coordinates());
    Matrix M2(MPoints);

    Matrix MT = M2.transpose();
    const vector<Point> lambda = MT.get_kernel();
    
    for(int i = 0; i < 3; ++i) {
        vector<mpq_class> pointLambda = lambda[0].get_coordinates();
        Point p = Point(B[i],true).multiply(pointLambda[i]);
        B[i]=p.get_coordinates();
    }
    
    Matrix Bt= Matrix(B).transpose();
    
    
    //vertices of polytope = kernel Bt, if necessary change row for all ones
    Matrix At = Matrix(Bt.get_kernel());
    Matrix A = At.transpose();
    
    bool ones = false;
    int ones_index = 0;
    for(int i = 0; i < 5; ++i) {
        bool ok = true;
        for(int j = 1; j < 8 and ok; ++j) {
            if (A.get(i,j) != A.get(i,j-1)) ok = false;
        }
        if (ok) {
            ones = true;
            ones_index = i;
        }
    }
    vector<vector<mpq_class>> matA = A.getMatrix();
    if (ones) {
        swap(matA[ones_index],matA[4]);
    }
    for(int i = 0; i < 8; ++i) {
        matA[4][i] = 1;
    }
    A = Matrix(matA);
    
    return A;
}



