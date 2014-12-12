'''
This routine compute all the matrix which have the same columns as A, up to permutations.
    @param: A matrix A.
    @return: A list of matrix containing all the possible columns-permuted matrices of A.'''
from numpy import zeros
from itertools import permutations
def perColumns(M):   
    (n,m) = M.shape    
    Sm=list(permutations(xrange(m)))
    all_permutations = []
    for alpha in Sm:
        A = zeros((n,m))
        for i in xrange(0,m):
            A[0:n,i] = M[0:n,alpha[i]]
        all_permutations.append(A)
    return all_permutations
