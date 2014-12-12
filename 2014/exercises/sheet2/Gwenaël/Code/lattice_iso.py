'''
This routine check if for any two given matrices P and Q with integer coefficients, it exists a matrix A with integer coefficient such that P = AQ' and det(A) = 1, 
where Q' is a matrix with the same columns as Q, but possibily permuted.
    @param: Two matrices P and Q of arbitrary dimension. 
    @return: a matrix A which satisfies the above properties (if it exists) ; 0 otherwise.'''
    
from numpy import *
from comatrix import *
from perColumns import *
from numpy.linalg import inv

def lattice_iso(P,Q):  
    equiv = 0
    (n,m) = Q.shape
    all_permutations = perColumns(Q)
    
    if abs(Determinant(P)) != abs(Determinant(Q)):
        print('The two matrices cannot be lattice isomorphic because |det(P)| != |det(Q)|') 
        return 0
    
    for Qp in all_permutations:
        coQ = comatrix(Qp)
        det = Determinant(Qp)
        A = dot(P,transpose(coQ))
        equiv = 1
        for i in xrange(0,n):
            if equiv == 0: break
            for j in xrange(0,m):
                if A[i,j] % det != 0:
                    equiv = 0
                    break
        if equiv == 1:
            print('The two matrices are lattice isomorphic.')
            A = A*(1.0/det)
            return(A)
    
    print('The two matrices are not lattice isomorphic.')
    return 0

