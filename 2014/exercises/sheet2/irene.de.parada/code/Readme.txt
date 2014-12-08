Function test_chi is given the name of a file containing in each row three file names.
 
First one corresponds to the file in which the homogeneous coordinates of the points
 of the first tetrahedron are stored in rows. Writing these coordinates in columns 
 we define matrix P.

Second one corresponds to the file in which the homogeneous coordinates of the points
 of the second tetrahedron are stored in rows. Writing these coordinates in columns 
 we define matrix Q.

The third file name corresponds to the output file in which the matrices defining a
 lattice isomorphism are written.

When reading the matrices it is checked that all entries are integers that they define
 a 4x4 matrix. It is also checked that any column is the zero vector and the homogeneous
 coordinates are transformed so that the first entry is a 1. Exceptions are raised for
 these cases.

Then we compute A= Q*P^{-1} for all matrices P with a permutation  of the columns of P.
 If A is a matrix of integers defining a lattice isomorphism and A has not been stored 
 as a solution then we add it to solutions.

** Instructions
Run the .script from the command line using 'sudo sage TestChi.sage'. Super user mode is 
needed to be able to write to the output files.

Upon start, the program will read the required tests and file-names from a file called 
'test_chi_input.txt' which must be located on the same folder.