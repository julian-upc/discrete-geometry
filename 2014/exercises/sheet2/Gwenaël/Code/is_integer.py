'''
This routine check if a given element is a integer.
    @param: an int, a long, a float, a complex or an array A.
    @return: True if A represent an integer number or is an array which contains only integer numbers ; False otherwise.
'''
from types import*
def is_integer(A):
    if(type(A) in (int, long)): return True
    if(type(A) in (float, complex)): return A == int(A)
    (n,m) = A.shape
    for i in xrange(0,n):
        for j in xrange(0,m):
            if A[i,j] != int(A[i,j]): return False
    return True
            
