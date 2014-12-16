'''
This routine compute the determinant of a given matrix M.
    @param: A matrix A.
    @return: The cofactor matrix of A.'''
from numpy import *
def Determinant(M):
    (n,m) = M.shape
    if((n,m)==(1,1)): 
        return(M[0,0])
    else:
        det = 0
        for j in xrange(0,n):
            A = delete(M,(0), axis=0)
            A = delete(A,(j), axis=1)
            det = det + M[0,j]*Determinant(A)*((-1)**(j+2))
        return det
