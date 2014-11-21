#!/usr/bin/python
from numpy import linalg
import numpy as np
from itertools import permutations 
from tetrahedron import calcA, determinant, cofactor, divisibility, divide, permuta, compara

def main(matrix_P,matrix_Q,matrix_A):#some main to test the function
	if testchi(matrix_P,matrix_Q,matrix_A)==True:
		print('Yes!')
	else:
		print('NO!')

def testchi (matrix_P,matrix_Q,matrix_A):#it does no checkings, for personal use only!
	P=np.loadtxt(matrix_P)
	Q=np.loadtxt(matrix_Q)
	A=np.loadtxt(matrix_A)
	P=np.transpose(P)
	Q=np.transpose(Q)
	d=determinant(P)
	alA=calcA(P,Q,d)
	if compara(A,alA)== True:
		return True
	else:
		return False

