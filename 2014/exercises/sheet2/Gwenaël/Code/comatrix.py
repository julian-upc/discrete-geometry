'''
This routine compute the cofactor matrix of a given matrix M.
    @param: A matrix A.
    @return: The determinant of A.'''
from numpy import *
from Determinant import *
def comatrix(M):
    (n,m) = M.shape
    CoM = ones((n,m), dtype=int16)
    for i in xrange(0,n):
        for j in xrange(0,m):
            A = delete(M,(i), axis=0)
            A = delete(A,(j), axis=1)
            CoM[i,j] = Determinant(A)*(-1)**(j+i+2)
    return CoM

