//
//  Matrix.cpp
//  PolytopeProject
//


#include "Matrix.hpp"

Matrix::Matrix(vector<Point> points){
    for(int i=0; i<points[0].get_coordinates().size(); i++){
        matrix.push_back(vector<mpq_class>());
    }
    
    for(int i=0; i<points.size(); i++){
        vector<mpq_class> coord = points[i].get_coordinates();
        for(int j=0; j<coord.size(); j++){
            coord[j].canonicalize();
            matrix[j].push_back(coord[j]);
        }
    }
}

Matrix::Matrix(vector<mpq_class> &array){
    for(int i=0; i<array.size(); i++){
        array[i].canonicalize();
        vector<mpq_class> vec = {array[i]};
        matrix.push_back(vec);
    }
}

Matrix::Matrix(vector<vector<mpq_class>> &vec): matrix(vec){}

void Matrix::multiplyRow(vector<mpq_class> &row,mpq_class num){
    for(int i=0; i<row.size(); i++){
        mpq_class multiplied = row[i]*num;
        multiplied.canonicalize();
        row[i] = multiplied;
    }
}
void Matrix::subtractRows(const vector<mpq_class> high_row,vector<mpq_class> &low_row, mpq_class scalar){
    for(int i=0; i<high_row.size(); i++){
        mpq_class subtrahend = high_row[i]*scalar;
        subtrahend.canonicalize();
        mpq_class difference = low_row[i]-subtrahend;
        difference.canonicalize();
        low_row[i]=difference;
    }
}


vector<Point> Matrix::get_kernel(){
    vector<Point> base_kernel;
    
    rref();
    
    for(int i = rowDim(); i < colDim(); i++){
        //Makes a vector of size colDim, with all zeroes
        vector<mpq_class> vec(colDim(),0);
        
        //Copy the elements from column i+1 to the vector
        for(int j=0; j<rowDim(); j++){
            vec[j] = matrix[j][i];
        }
        
        //Add the -Identity matrix
        vec[i] = -1;
        
        Point point(vec,true);
        base_kernel.push_back(point);
    }
    
    return base_kernel;
}

mpq_class Matrix::get(int row,int col) const { return matrix[row][col];}
int Matrix::rowDim() const{return (int)matrix.size();}
int Matrix::colDim() const{return (int)matrix[0].size();}
vector<vector<mpq_class>> Matrix::getMatrix() const {return matrix;}

//Multiplies this matrix with B and returns the result
Matrix Matrix::multiply(Matrix B) const {
    if(colDim()!=B.rowDim()){ throw runtime_error("Wrong Matrix Dimension");}
    
    vector<vector<mpq_class>> C;
    
    for(int i=0; i<rowDim(); i++){
        vector<mpq_class> row;
        for(int j=0; j<B.colDim(); j++){
            mpq_class sum=0;
            for(int k=0; k<colDim();k++){
                sum = sum + matrix[i][k]*B.get(k,j);
            }
            sum.canonicalize();
            row.push_back(sum);
        }
        C.push_back(row);
    }
    return Matrix(C);
}

void Matrix::print() const{
    cout<<endl;
    for(int i=0; i<matrix.size(); i++){
        for(int j=0; j<matrix[i].size(); j++){
            cout << matrix[i][j].get_d() << " ";
        }
        cout<<endl;
    }
}

Point Matrix::pivot(int column_index){
    mpq_class max=abs(matrix[column_index][column_index]);
    int max_row = column_index;
    for(int i=column_index+1; i<rowDim(); i++){
        mpq_class next_num = abs(matrix[i][column_index]);
        if(next_num > max){
            max = next_num;
            max_row=i;
        }
    }
    if(max==0) throw runtime_error("Matrix is singular!");
    
    mpq_class pivotValue = matrix[max_row][column_index];
    multiplyRow(matrix[max_row], 1/pivotValue);
    
    bool swapped = false;
    if(max_row != column_index){
        swapped = true;
        swap(matrix[column_index],matrix[max_row]);
    }
    
    
    return Point(vector<mpq_class> {pivotValue},swapped);
}

//Note: I also calculate the lower triangle (all zeroes), this can be optimised if necessary
mpq_class Matrix::rref(){
    mpq_class determinant = 1;
    
    int noSwaps=0;
    
    for(int k=0; k<min(rowDim(),colDim()); k++){
        //Find pivot
        Point pivotPoint = pivot(k);
        
        determinant*=pivotPoint.get_coordinates()[0];
        determinant.canonicalize();
        
        noSwaps += pivotPoint.get_sign();
        
        //Do for all the rows below the pivot
        for(int i=0; i<rowDim(); i++){
            if(i!=k && matrix[k][k]!=0){
                mpq_class multiplier = matrix[i][k]/matrix[k][k];
                multiplier.canonicalize();
                
                subtractRows(matrix[k], matrix[i], multiplier);
            }
        }
    }
    
    if(noSwaps%2!=0) determinant*=-1;
    
    return determinant;
}

mpq_class Matrix::determinant(){
    if(rowDim()==2 && colDim()==2){
        mpq_class pos = matrix[0][0]*matrix[1][1];
        mpq_class neg = matrix[0][1]*matrix[1][0];
        pos.canonicalize();
        neg.canonicalize();
        mpq_class result = pos-neg;
        return result;
    }
    return(rref());
}

Matrix Matrix::transpose() const{
    vector<vector<mpq_class>> nextMatrix;
    for(int i=0; i<colDim(); i++){
        vector<mpq_class> row;
        for(int j=0; j<rowDim(); j++){
            row.push_back(matrix[j][i]);
        }
        nextMatrix.push_back(row);
    }
    return Matrix(nextMatrix);
}

void Matrix::testMatrix(){
    Point p1(vector<double>{1,2},true);
    Point p2(vector<double>{3,7},true);
    Point p3(vector<double>{5,1},true);
    
    Matrix A(vector<Point> {p1,p2,p3});
    
    Point q1(vector<double>{1,3,5},true);
    Point q2(vector<double>{2,4,6},true);
    Point q3(vector<double>{7,3,1},true);
    
    Matrix B(vector<Point> {q1,q2,q3});
    cout << "determinant" << B.rref().get_d() << endl;
    
    Point r1(vector<double>{0,2},true);
    Point r2(vector<double>{1,7},true);
    cout << Matrix(vector<Point>{r1,r2}).determinant().get_d() << endl;
    
    vector<mpq_class> vec3{mpq_class(19,3102),mpq_class(128,127893)};
    Matrix D(vec3);
    
    A.print();
    B.print();
    Matrix C=A.multiply(B);
    C.print();
    D.print();
    
    //C.multiply(D).print();

    print();
    mpq_class t = mpq_class(13,17);
    multiplyRow(matrix[0],t);
    print();
    subtractRows(matrix[0],matrix[1],matrix[1][0]/matrix[0][0]);
    print();
    
    Point a1(vector<mpq_class>{mpq_class(1,2),mpq_class(1,3)},true);
    Point a2(vector<mpq_class>{mpq_class(1,3),mpq_class(1,5)},true);
    Point a3(vector<mpq_class>{mpq_class(19,17),mpq_class(13,12)},true);
    Point a4(vector<mpq_class>{mpq_class(111,3334),mpq_class(1,120137280713027)},true);
    Matrix M(vector<Point>{a1,a2,a3,a4});
    
    M.print();
    //M.pivot(0);
    //M.print();
    vector<Point> vec = M.get_kernel();
    Matrix N = Matrix(vec);
    N.print();
    M.multiply(N).print();
    N.print();
    Matrix O = N.transpose();
    O.print();
    O.transpose().print();
    
}
