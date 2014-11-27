In this file you will find two files :
- chi.py
- test_chi.py

As the name says it the function test_chi will test the function chi.
It takes in argument a txt file where the lines have the following form : 
	matrix_P1.txt matrix_P2.txt result_P1_P2.txt
Each line contains a set of argument for the function chi

The function chi takes three argument 
	- matrix_P and matrix_Q are the txt file containing the matrices P and Q
	- matrix_A is where the function is going to write the matrix A it founds

The matrices P and Q are written in the following way in the txt file : 
	1 0 0 0 
	1 1 0 0 
	1 0 1 0 
	1 2 3 5

And the chi function make some test on the entries before calling the matrix computation function : 
	- the matrices must be big enough, 4x4 or more and it will take the upper left 4x4 matrix
	- the matrices must be integer matrices
	- the matrices must represent a tetrahedron, so be of full rank.
	- the matrices must have the same determinants in absolute value, otherwise they can not be lattice equivalent. 

