# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 20:30:45 2014

@author: Cassandra de Cock de Rameyen
"""

import numpy as np
import numpy.linalg as lin
from itertools import permutations

def chi(matrix_P, matrix_Q, matrix_A):

    P = np.transpose(np.loadtxt(matrix_P))
    Q = np.transpose(np.loadtxt(matrix_Q))
    
    #flags for testing the entries
    test_size = True
    test_entry = True

    # testing the size of the entries, 
    # if bigger than 4x4 take the 4x4 upper left part of the matrix
    # if to small return that the test is false
    if (P.shape[0] >= 4 and P.shape[1] >= 4):
        P = P[0:4,0:4]
    else:
        test_size = False

    if (Q.shape[0] >= 4 and Q.shape[1] >= 4):
        Q = Q[0:4,0:4]
    else:
        test_size = False 

    # test if all the entries are integer,
    # if Q and P have the same determinant in absolute value
    # if P an Q are of full rank
    if test_size :
        for i in range(4):
            for j in range(4):
                if P[i,j] != int(P[i,j]):
                    test_entry = False
                if Q[i,j] != int(Q[i,j]):
                    test_entry = False
        if abs(determinant(P,4)) != abs(determinant(Q,4)):
            test_entry = False
        
        if lin.matrix_rank(P)!=4 or lin.matrix_rank(Q)!=4 :
            test_entry = False

    # Compute the matrix A if P and Q passed all the test
    if test_size & test_entry :
        file = open(matrix_A,'w')
        list_A = matrix_computation(P, Q, file)
        file.close()
    else:
        print('Wrong Entries')

# Computation of matrix A
def matrix_computation(P, Q, file):

	inttest = True
	list_A = []
      
      # for each permutation of P, compute it's adjugate
      # check after the product Q*adj(P) if the entries are divisible by det(P)
      # check if the matrix was already known
      # if it didn't find any matrix it return the zero matrix
	for p in permutations(range(4)):
     
		Pnew = np.transpose(np.array([P[:,p[0]], P[:,p[1]], P[:,p[2]], P[:,p[3]]]))		
		A = np.dot(Q,adjugate(Pnew))
		det = determinant(Pnew,4)

		for i in range(4):
			for j in range(4):
				if abs(A[i,j] % det) == 0:
					A[i,j] = A[i,j]/det
				else:
					inttest = False
     
				if abs(A[i,j]) == 0:
					A[i,j] = 0

		if inttest and not(any((A == X).all() for X in list_A)):
			file.write('A = \n')
			np.savetxt(file,A, fmt='%s')
			list_A.append(A)

		inttest= True

	if len(list_A)== 0:
		A = np.zeros((4,4))
		file.write('A = \n')
		np.savetxt(file,A, fmt='%s') 


def adjugate(P):
	cofactor = np.zeros((4,4))
	for i in range(4):
		for j in range(4):
			minor = np.delete(np.delete(P,i,0),j,1)
			cofactor[i,j] = (-1)**(i+j)*determinant(minor,3)
	return np.transpose(cofactor)

def determinant(M,n):
	if n == 1:
		det = M[0,0]
	if n == 2:
		det = M[0,0]*M[1,1]-M[0,1]*M[1,0]
	else:
		det = 0
		for i in range(n):
			minor = np.delete(np.delete(M,i,0),0,1)
			det += (-1)**(i)*determinant(minor,n-1)*M[i,0]
	return det 
