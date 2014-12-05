# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

# Given a letter X \in {A,D,E,F,H}, a positive integer n and as optional input a vector v
# with n+1 coordinates (the first must be zero), the function
#
# v_orbit_cardinality (x, n [, v])
#
# returns the cardinality of the complete set of the reflections of v under X.
#
#If v is not specified, the function choose a vector not on any reflecting hyperplane.
# -----------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------

import numpy as np
import random
import time


def matrix_A (n):
	n = n + 1
	matrix_adj = np.zeros([n-1,n+1])
	for i in range(n-1):
		matrix_adj[i][i+1] = 1
		matrix_adj[i][i+2] = -1
	return matrix_adj


def matrix_D (n):
	matrix_adj = np.zeros([n,n+1])
	for i in range(n-1):
		matrix_adj[i][i+1] = 1
		matrix_adj[i][i+2] = -1
	matrix_adj[n-1][n-1] = 1
	matrix_adj[n-1][n] = 1
	return matrix_adj

def matrix_E (n):
	matrix_adj = np.zeros([n,n+1])
	for i in range(n-2):
		matrix_adj[i][i+1] = 1
		matrix_adj[i][i+2] = -1
	matrix_adj[n-2][n-2] = 1
	matrix_adj[n-2][n-1] = 1
	for i in range(n):
		matrix_adj[n-1][i+1] = 1
	return matrix_adj


def ref (v, a): #reflection of v for the hyperplane orthogonal to a (going through the origin)
	ref_a_v=np.zeros([len(a)])
	for i in range(len(a)):
		ref_a_v[i] = v[i] - (2 * a[i] * np.dot(v,a) / np.dot(a,a) )
	return ref_a_v

def reflected_vectors(X):
	reflected_vectors = []
	reflected_vectors.append(X[0])
	for vector in reflected_vectors:
		for hyperplane in X:
			vector_reflection = ref(vector, hyperplane)
			k = 0
			for vector_2 in	reflected_vectors:
				if equal_vectors(vector_2, vector_reflection) == True:
					k += 1
					break
			if k == 0:		
				reflected_vectors.append(vector_reflection)	
	return reflected_vectors

def vector_not_in_hyperplanes(X):
	v = np.zeros(len(X[0]))
	v[0] = 0
	Y = reflected_vectors(X)
	for i in range(len(Y[0])-1):
		v[i+1] = random.randint(1,100)
	for j in range(len(Y)):
		if abs(np.dot(Y[j],v)) < epsilon:
			vector_not_in_hyperplanes(X)
	return v

def equal_vectors(X,Y): #returns true if the vectors X and Y are equals
	if len(X) != len(Y):
		return False
	else:
		for i in range(len(X)):
			if abs(X[i]-Y[i]) > epsilon:
				return False
	return True


# #Cardinality of the Complete set of reflecting hyperplanes
def v_orbit_cardinality (x, n, v=False):
	if x == "A":
		X_n = matrix_A(n)
	
	elif x == "D":
		X_n = matrix_D(n)
	
	elif x == "E" and n == 8:	
		X_n = matrix_E(n)

	elif x == "F" and n == 4:
		X_n = [[0,1,-1,0,0], [0,0,1,-1,0], [0,0,0,1,0], [0,-1,-1,-1,-1]]
	
	elif x == "H" and n == 3:
		X_n = [[0,2,0,0], [0, -3-sqrt(5), 2, -1-sqrt(5)], [0,0,0,2]]
	
	elif x == "H" and n == 4:
		X_n = [[0,-3-sqrt(5),2,2,2], [0,-1,1,0,0], [0,0,-1,1,0], [0,0,0,-1,1]]
	else:
		print "Select x\in {A,D,E,F,H} and if \n x = E then n = 8\n x = F then n = 4 \n x = H then n = 3 or n=4"
		return False

	if v == False:
			v = vector_not_in_hyperplanes(X_n)

	print X_n
	# print "vector= ", v
	vectors = []
	vectors.append(v)
	for vector in vectors:
		for hyperplane in X_n:
			vector_reflection = ref(vector, hyperplane)
			k = 0
			for vector_2 in vectors:
				if equal_vectors(vector_2, vector_reflection) == True:
					k += 1
			if k == 0:
				vectors.append(vector_reflection)
	# print return
	return len(vectors)



#------------------------------------------------------------------------------------------------------

epsilon = 0.00001
start_time = time.clock()

# print "Cardinality = ", v_orbit_cardinality("D", 4, [0,64,72,63,31])
print "Cardinality = ", v_orbit_cardinality("A", 6)

print "Running time = ", time.clock() - start_time, "seconds"


# ----------------------------------------------------------  TO IMPLEMENT  ----------------------------------- 
# Matrices H (with sqrt(5))

# ----------------------------------------------------------  TESTS ----------------------------------- 


# Order of An = (n+1)!
	# if n=3 order of An= 4! = 24 (checked) T= 0,015 seconds
	# if n=4 order of An= 5! = 120 (checked) T= 0,21 seconds
	# if n=5 order of An= 6! = 720 (checked) T= 2,09 seconds
	# if n=6 order of An= 7! = 5040 (checked) T = 417,70 seconds
	# if n=7 order of An= 8! = 40320


# Order of Dn = 2^(n-1) n!
	# if n=3 order of Dn = 2^2 * 3! = 24 (checked)
	# if n=4 order of Dn = 2^3 * 4! = 192 (checked) 0.44 sec
	# if n=5 order of Dn = 2^4 * 5! = 1920 (checked) 51 sec
	# if n=6 order of Dn = 2^5 * 6! = 23040

# Order of F4 = 1152 (checked)

