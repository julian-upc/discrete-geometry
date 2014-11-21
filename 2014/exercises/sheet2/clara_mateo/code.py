import numpy as np
import csv
from itertools import permutations

# Auxiliar Functions -------------------------------------------------------

def matrix_multiplication(X,Y): #it multiplies two matrices
	result = np.dot(X,Y)
	return result

def matrix_transpose(X): #it trasposes a matrix
	matrix_trans = np.zeros([len(X),len(X[0])])
	for i in range(len(X)):
		for j in range(len(X[0])):
			matrix_trans[i][j] = (X[j][i])
	return matrix_trans

def matrix_determinant(X):
	return np.linalg.det(X)

def matrix_inverse(X):
	Y = np.linalg.inv(X)
	return Y

def matrix_zeros(X): #Check if there is any "zero" in a matrix
	for j in range(len(X)):
		for i in range(len(X[0])):
			if abs(X[j][i]) < 0.000001:
				X[j][i] = 0.
	return X

def matrix_integer(X): #it checks if the matrix has int entries
	for j in range(len (X)):
		for i in range(len (X[0])):
			a = X[j][i]
			if (a).is_integer() == False:
				return False
	return True

def matrix_cofactors(X): # Compute the matrix of cofactors
	if len(X) ==  len(X[0]):
		matrix_adj = np.zeros([len(X),len(X[0])])
		for j in range(len(X)): #rows
			for i in range(len(X[0])): #columns
				matrix_aux = np.zeros([len(X)-1,len(X[0])-1])
				matrix_aux = np.delete(np.delete(X,j,0),i,1) #j delete row, i delete column
				matrix_adj[j][i] = pow(-1,i+j) * matrix_determinant(matrix_aux)
				# print "j= ", j, "i= ", i, X[j][i], "matrix_aux= ", matrix_aux, "det", matrix_determinant(matrix_aux)
		return matrix_adj
	else:
		print "The matrix", X, "is not square!!"
		return X

def matrix_adjugate(X): #Compute the adjugate matrix
	matrix_adj = matrix_cofactors(X)
	return matrix_transpose(matrix_adj)

def equal_matrices(X,Y): #Check if two matrices are equal
	if len(X) == len(Y) and len(X[0]) == len(Y[0]):
		for j in range(len(X)):
			for i in range(len(X[0])):
				if abs(X[j][i] - Y[j][i]) > 0.000001:
					return False
		return True
	else:
		return False
		
def matrix_id(X): #Give the Id matrix of the sime sice as X (if this one is square)
	if len(X) ==  len(X[0]):
		matrix_id = np.zeros([len(X),len(X[0])])
		for j in range(len(X)):
			matrix_id[j][j] = 1.
		return matrix_id
	else:
		print "is not square", X
		return X 

# Exercise Function itself ------------------------------------------------------------------

def homogeneous_map(X):
	# Files that we want to read (with the matrices P and Q)
	files = [X[0], X[1]]
	zero_def = 1e-6 #A value less than this will be considered as 0.0

	#Reading and saving the matrices...
	i = 0
	for file in files:
		f = open(file, 'r')
		matrix = []
		matrix = [ map(int,line.split(' ')) for line in f ]
		if i == 0:
			matrix_1 = matrix_transpose (matrix)
		else:
			matrix_2 = matrix_transpose(matrix)
		i += 1

	# det P = \pm det Q, otherwise there's no affine transformation
	if abs(abs(matrix_determinant(matrix_1)) - abs(matrix_determinant(matrix_2))) < zero_def:


		if matrix_integer(matrix_1) == False or matrix_integer(matrix_2) == False:
			print "The initial matrices must have integer coefficients"
			np.savetxt(X[2], np.zeros([len(matrix_1),len(matrix_1[0])]), fmt="%i")
		else:
			set_of_matrices = [] #This will be where the final result will be saved

			#consider all the permutations of matrix_1
			for p in permutations(range(len(matrix_1))):
				matrix_1_perm_aux = np.zeros([len(matrix_1),len(matrix_1[0])])
				matrix_1_perm = np.zeros([len(matrix_1),len(matrix_1[0])])
				for i in range(len(matrix_1)): #this will give the permuted matrix
					matrix_1_perm_aux[p[i]] = matrix_transpose(matrix_1)[i] 
				matrix_1_perm = matrix_transpose(matrix_1_perm_aux)

				#We compute adjugate(P)*Q
				matrix_adj = matrix_adjugate(matrix_1_perm)
				matrix_X = matrix_multiplication(matrix_adj, matrix_2)
		 		
		 		#Now, we check if all the entries of adjugate(P)*Q (matrix_X) are divisible by det(matrix_1_perm)
		 		num_of_no_int = 0
		 		for j in range(len(matrix_X)):
		 			for i in range(len(matrix_X[0])):
		 				det_matrix = matrix_determinant(matrix_1_perm)
		 				mod = abs(matrix_X[j][i]) % abs(det_matrix)
		 				if det_matrix != 1:
			 				if abs(matrix_X[j][i]) > 0.000001: 
			 					if abs(matrix_X[j][i]) < (abs(det_matrix)-0.5):
			 						num_of_no_int += 1
			 						# print "0. j= ", j, "i= ", i, "Xji= ", matrix_X[j][i], "det= ", det_matrix
			 					else:
			 						if mod > 0.5 and mod < (abs(det_matrix)-0.5):
				 						num_of_no_int += 1
				 						# print "1. j= ", j, "i= ", i, "Xji= ", matrix_X[j][i], "det= ", det_matrix, "mod= ", mod, (abs(matrix_X[j][i]) % abs(det_matrix)), num_of_no_int
		 				

		 		if num_of_no_int == 0: #if 1/det(P) adjugate(P)*Q \in M(\Z)
		 			matrix_trans = matrix_zeros(matrix_X / det_matrix)

					if len(set_of_matrices) == 0:
						set_of_matrices.append(matrix_trans)
					else: #Check if there are no equal matrices or the identity matrix
						k = 0
						for num in range(len(set_of_matrices)):
							if equal_matrices(matrix_trans, matrix_id(matrix_trans)) == True:
								k += 1
							if equal_matrices(set_of_matrices[num],matrix_trans) == True:
								k += 1
						if k == 0:
							set_of_matrices.append(matrix_trans)


			#Saving it in a new file
			if len(set_of_matrices) == 0:
				np.savetxt(X[2], np.zeros([len(matrix_1),len(matrix_1[0])]), fmt="%i")
			else:
				np.savetxt(X[2], set_of_matrices[0], fmt="%f",  newline='\n')
				for i in range(1, len(set_of_matrices)):
					with open(X[2],'a') as f_handle:
		   	 			f_handle.write("\n")
		   	 			np.savetxt(f_handle, set_of_matrices[i], fmt="%f", newline='\n')
			
	else:
		np.savetxt(X[2], np.zeros([len(matrix_1),len(matrix_1[0])]), fmt="%i")
		
		return 0

# Test-Chi Function ------------------------------------------------------------------

def test_chi(X):
	matrix =[]
	with open(X,'r') as f:
	    for line in f:
	        for word in line.split():    
	           matrix.append(word)
	a = len(matrix)
	while (a > 2):
		matrices_files = [matrix[0], matrix[1], matrix[2]]
		del matrix[0]
		del matrix[0]
		del matrix[0]
		a = a - 3
		homogeneous_map(matrices_files)
	return 0


print test_chi('prova.txt')
