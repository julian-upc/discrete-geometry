#!/usr/bin/python
#-*- coding: utf-8-*-
from numpy import linalg
import numpy as np
from itertools import permutations 

def main(matrix_P, matrix_Q, matrix_A):#loads matrices and checks for bad inputs (bigger/smaller matrices) and zero det matrices
	P=np.loadtxt(matrix_P)
	Q=np.loadtxt(matrix_Q)
	P=np.transpose(P)
	Q=np.transpose(Q)
	f = open(matrix_A,'w')
	if P.shape!=(4,4) or Q.shape!=(4,4):
		#error message
		f.write('matrix P or Q do not have the appropriate dimensions')
	else:
		if (entera(P)==False or entera(Q)==False):
			f.write('matrix P or Q do not have integer values')
		else:
			d=determinant(P)
			dq=determinant(Q)
			if (d==0 or dq==0):
				#error message
				f.write('P or Q have 4 coplanar points')
			else:
				alA=calcA(P,Q,d)
				writeA(alA, f)
	f.close()

def entera(P):#checks if the entry has integer values
	for i in range(4):
		for j in range(4):
			if P[i,j]!= int(P[i,j]):
				return False
	return True

def writeA(alA, f):#writes the matrices found in file f
	if not alA:
		m=np.zeros((4,4))
		np.savetxt(f,m,fmt='%s')
	else:
		for a in alA:
			f.write('A=\n')
			np.savetxt(f,a,fmt='%s')
			f.write('\n')
	

def calcA(P,Q,d):#given two matrices and determinant of the first computes possible transformations
	perQ=permuta(Q)
	alA=[]
	P2=cofactor(P)
	for q in perQ:
		A= np.dot(q,P2)
		if divisibility(A,d)==True:
			A=divide(A,d)
			if compara(A,alA)==False:
				alA.append(A)
	return alA


def determinant(P):#computes determinant of a matrix
	if P.shape[0]==2:
		det=P[0,0]*P[1,1]-P[1,0]*P[0,1]
	else:
		det=0
		for i in range(P.shape[0]):
			P1=np.delete(P,i,0)
			P1=np.delete(P1,0,1)
			p1=determinant(P1)
			det=det+((-1)**i)*P[i,0]*p1
	return det

def cofactor(P):#compute matrix of cofactors transposed
	P2=np.zeros((4,4))
	for i in range(4):
		for j in range(4):
			P1=np.delete(P,i,0)
			P1=np.delete(P1,j,1)
			P2[i,j]=((-1)**(i+j))*determinant(P1)

	return np.transpose(P2)

def divisibility(A,d):#checks if matrix is divisible by integer
	for i in range(4):
		for j in range(4):
			if abs(A[i,j]%d)!=0:
				return False
	return True

def divide(A,d):#divides matrix by integer
	for i in range(4):
		for j in range(4):
				A[i,j]/=d
	return A
 
def permuta(Q):#returns a list with all the possible column permutations of the matrix Q
	P=[0,1,2,3]
	per=list(permutations(P)) 
	perQ=[]
	for p in per:
		Q1=np.zeros((4,4))
		for i in range(4):
			for j in range(4):
				Q1[i,j]=Q[i,p[j]]
		perQ.append(Q1)
	return perQ

def compara(A, alA):#checks if a matrix is in a list of matrices
	for a in alA:
		comp= True
		for i in range(4):
			for j in range(4):
				if a[i,j]!=A[i,j]:
					comp=False
		if comp== True:
			return True
	return False