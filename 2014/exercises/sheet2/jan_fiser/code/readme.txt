TASK

The program test_chi decides whether two 3-dimensional lattice tetrahedra P and Q are lattice equivalent, i. e., whether there exists an affine map f(X) such that f(P) = Q.

------------------------------------------------------------

REPRESENTATION OF THE ALGEBRAICAL OBJECTS

The tetrahedron P (resp. Q) is represented by a 4x4 matrix (whose columns are vertices of P (resp. Q) in homogeneous coordinates). 

The affine map f(x) = A'x + b is represented by a 4x4 matrix A:

| 1    0   0   0 |
| b_1            | := A
| b_2      A'    |
| b_3            |

where A' is 3x3 matrix with determinant +-1.

------------------------------------------------------------
  
INPUT & OUTPUT

INPUT: Tetrahedra P,Q written in textfiles with names chosen by the user. Every file contains one matrix representing one polyhedron as follows:

p(1)_1 p(1)_2 p(1)_3 p(1)_4
p(2)_1 p(2)_2 p(2)_3 p(2)_4
p(3)_1 p(3)_2 p(3)_3 p(3)_4
p(4)_1 p(4)_2 p(4)_3 p(4)_4

where row i represents i-th vertex of the polyhedron.
 

OUTPUT: If P is lattice eqiuvalent to Q, then test_chi retruns set of all affine maps corresponding to the equivalence, i. e., all matrices A such that AP = Q. Otherwise the program returns 0. The output is a textfile with a name chosen by the user.

------------------------------------------------------------

PROGRAM

0. Read matrices P,Q from the input files.

1. Check if they represent two tetrahydra:
	- P,Q are of the shape 4x4,
	- P,Q have integer entries,
	- P,Q are full-dimensional,
	- P,Q are invertible.

2. If det(P) <> +-det(Q) go to step 4.  #P and Q cannot be lattice equivalent

3. For every permutation of columns of P:
	if A=PQ^{-1} has integer entries, add A to the list of all solutions

4. Print all solutions (uniquely) to the output textfile.


------------------------------------------------------------

USED TECHNIQUES

Determinant is computed from a definition to keep it as integer value.

Inverse matrix is represented as a pair (adjungate matrix, determinant),
and then the multiplying A=PQ^{-1} is 
	- multiply P and adungate Q =:A*
	- if all entries A* are divisible by the determinant of Q, A=A*/det Q.
	  Otherwise this case has no solution.

------------------------------------------------------------

USER MANUAL

To the main directory put your input textfiles, e. g. matrix_P1.txt and matrix_P2.txt.

In the code on the lines 6, 7, 8 fill
	- name of the input file containing P (matrix_P1.txt),
	- name of the file input containing Q (matrix_P2.txt),
	- name of output file (e.g. result_P1_P2.txt).

Run test-chi.py.
	
