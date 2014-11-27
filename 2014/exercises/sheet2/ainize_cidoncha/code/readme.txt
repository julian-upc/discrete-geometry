 
   Function mylatt2.m
   -----------------

  The input for this function is two files containing matrices with the coordinates of a lattice tetahedron in homogeneous coordinates and a name for the new file is going to be created. The function computes the A matrices with integer coordinates that map the first matrix into the other one.


test_chi.m
   --------------------------------

  The arguments of this function are three files. The first two must contain a matrix with the coordinates of a lattice tetahedron in homogeneous coordinates and the third one a matrix that maps the first tetahedron into the second one.

The program runs mylatt2.m into the first two files and compares the result with the one entered in the third file. If any of the matrices of the output of mylatt2 coincides with the one in the result it returns 1. Otherwise, it returns a 0.