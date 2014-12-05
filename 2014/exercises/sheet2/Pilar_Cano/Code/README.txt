LatticeEquivalence is a program that compares two tetrahedrons 
in homogenous coordinates

---------------------------------------------------------------------
=LatEq=
The main function is called LatEq(P,Q,R), 
where P and Q are your matrices and R is a file where the 
program should print your answer.
There are 3 types of answer:
-If they arre Lattice equivalence is going to answer one or 
several non-zero matrices.
-If they are not Lattice Equivalence, it is going to answer the
zero-matrix.
-If one of your entries is not a tetrahedron, or it is not in the
lattice, is going to answer 'I only work with Lattice tetrahedron'.

-----------------------------------------------------------------------
=TestSuit=
For the test of the program we have the function TestSuit(P,Q,R)
where P and Q are files that my program is going to read 4x4 matrices
of the form
1  0  0  0
1  1  0  0 
1  0  1  0
1  2  3  5

And R is a file where my program is going to return the answer of my 
function LatEq, this function. 

--------------------------------------------------------------------------
*In this program we have other functions that would be of you interest:


=Perm=
The funtion Perm(A) works by giving him a 4x4 matrix and gives you all the
possible permutations of your columns.

---------------------------------------------------------------------------
=Det=
The function Det(A,n) computes the determinant of a nxn matrix, in the case 
when you are working with Lattices, python does not give you the exct
answer that should be an integer, well, this function does it!

---------------------------------------------------------------------------- 
=ComMa=
The function ComMa(A,T,n) compares if two nxn matrices are the same, it will
return 0 if not and 1 if they are the same.

---------------------------------------------------------------------------- 
=ElemSet=
The function ElemSet(A,S) checks if the nxn A matrix is an element of a 
list of nxn matrices S. 

----------------------------------------------------------------------------
=Div=
The function DIV(A,d), is a function that if you give it a nxn matrix and an
integer d, will answer if all the entries of A are divisible by d or not,
if the answer is positive will return 1, if not, 0.

-----------------------------------------------------------------------------
=Integ=
Ingte(A) is a function that tells you if a nxn matrix A have all integer 
entries.

------------------------------------------------------------------------------
**How to run it!

   This program is coded in Python, that is whay the extension of our code is
.py, so, if you want to calculate something you define your 'somethin' put the 
function(one of LatticeEquivalence has) with its respective entries at the end 
of LatticeEquivalence.py .
 
 Go to your terminal type cd and the directory where LatticeEquivalence.py is, 
enter and then type python LatticeEquivalence.py, enter.

HaveFun!

--------------------------------------------------------------------------------
If you have some problems runing it or some advises in order to improve it
please contact with:

Pilar Cano, at
pilukno@gmail.com
