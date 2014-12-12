'''
This routine test the function chi with some specifical tetrahedrons.
'''
from chi import *
import numpy as np

# Test 1 : as suggested in sheet 2.
chi("matrix_P1", "matrix_Q1","result_1")

# Test 2 : homemade example
Q = array([(12,5,4,1), (0,8,6,12), (3,3,0,1), (2,7,9,6)])
A = array([(1,-8,-3,4), (0,1,6,-2), (0,0,1,1), (0,0,0,1)]) # upper triangulare matrix of determinant 1.
P = dot(A,Q)
matrix_P2 = open("matrix_P2", 'w')
matrix_Q2 = open("matrix_Q2", 'w')
np.savetxt(matrix_P2,transpose(P),fmt='%s')
np.savetxt(matrix_Q2,transpose(Q),fmt='%s')
B = chi("matrix_P2", "matrix_Q2","result_2")
if A == B: print("Test successful")
else: print("Test failed")

# Test 2 : homemade example
Q = array([(12,5,4,1), (0,8,6,12), (3,3,0,1), (2,7,9,6)])
A = array([(1,-8,-3,4), (0,1,6,-2), (0,0,1,1), (0,0,0,1)])
P = dot(A,Q)
matrix_P2 = open("matrix_P2", 'w')
matrix_Q2 = open("matrix_Q2", 'w')
np.savetxt(matrix_P2,transpose(P),fmt='%s')
np.savetxt(matrix_Q2,transpose(Q),fmt='%s')
B = chi("matrix_P2", "matrix_Q2","result_2")
if A == B: print("Test successful")
else: print("Test failed")
chi("matrix_P3", "matrix_Q3","result_3")


