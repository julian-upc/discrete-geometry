'''
This routine take as an input two files which contains respectively the matrix representation, in the projective space of dimension 3, 
of two distinct lattices tetrahedrons P and Q. It write in a third file, passed as argument, the matrix representation of a lattice 
isomorphism A which send P to Q (if it exists). If such an isomorphism does not exist, the file is left empty.
    @param: P_array, Q_array and A_array, three string which are the pathname of the file containing the matricial representation of 
            A, P and Q (respectively)
    @output: none.
'''
import numpy as np
from Determinant import *
from is_integer import *
from lattice_iso import *

def chi(P_array, Q_array, output):
    P=np.loadtxt(P_array)
    Q=np.loadtxt(Q_array)
    P=np.transpose(P)
    Q=np.transpose(Q)
    file = open(output, 'w')
    
    assert P.shape == (4,4) and Q.shape == (4,4),"P and Q should be of dimension 4x4."
    assert is_integer(P) and is_integer(Q), "Either P or Q is not an matrix with integer coefficients." 
    assert Determinant(P) != 0 and Determinant(Q) != 0, "Either P or Q is a non-invertible matrix and does not represent a tetrahedron."
    
    A = lattice_iso(P,Q)
    file.write('A=\n')
    np.savetxt(file,transpose(A),fmt='%s')
    file.close()