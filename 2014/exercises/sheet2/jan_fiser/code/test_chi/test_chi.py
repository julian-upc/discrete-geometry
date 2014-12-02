import numpy as np

# Here write names of the files containing input matrices P and Q. 
#To the solution_file variable write a name of the file to which the solution will be written.

P = np.transpose(np.loadtxt('matrix_P5.txt'))
Q = np.transpose(np.loadtxt('matrix_P6.txt'))
solution_file = 'result_P5_P6.txt'


S4 = np.loadtxt('s4.txt') 
# s4.txt contains on its first 12 lines even permutations of S4 and on the remaining 12 lines odd permutations of S4
# this fact we use in computing the determinant of 4x4 matrix
S3 = np.loadtxt('s3.txt')
# s3.txt contains on its first 3 lines even permutations of S3 and on the remaining 3 lines odd permutations of S3
# this fact we use in computing the determinant of 3x3 matrix

def first_checking(M): 

	# M matrix
	# return 1 if the given matrix M is of the correct 4x4 shape, full-dimensional and invertible
	# return 0 otherwise

	if (M.shape != (4,4)):
		print('P or Q is not a 4x4 matrix.')
		return(0)
	for i in range(4):
		for j in range(4):
			if (M[i][j] != int(M[i][j])):
				print('P or Q contains non-integer entry.')
				return(0)
	if (determinant4(M) == 0):
		print('P or Q is not full-dimensional.')
		return(0)
	return(1)

def determinant4(M): 

	# M 4x4 integer matrix
	# return determinant of M
	
	d = 0
	for i in range(12):
		d += M[0][S4[i][0]-1]*M[1][S4[i][1]-1]*M[2][S4[i][2]-1]*M[3][S4[i][3]-1]	
	for i in range(12,24):
		d -= M[0][S4[i][0]-1]*M[1][S4[i][1]-1]*M[2][S4[i][2]-1]*M[3][S4[i][3]-1]	
	return(d)	

def determinant3(M): 

	# M 4x4 integer matrix
	# return determinant of M

	d = 0
	for i in range(3):
		d += M[0][S3[i][0]-1]*M[1][S3[i][1]-1]*M[2][S3[i][2]-1]	
	for i in range(3,6):
		d -= M[0][S3[i][0]-1]*M[1][S3[i][1]-1]*M[2][S3[i][2]-1]	
	return(d)	
	
def cofactor(M,x,y):

	# M 4x4 matrix
	# return (x,y)-cofactor of M (i.e., delete row x and column y from M)

	C = np.empty(shape=[3, 3])
	i = 0
	for ii in range(4):
		if (ii != x):
			j = 0
			for jj in range(4):
				if (jj != y):
					C[i][j] = M[ii,jj]
					j += 1
			i += 1		
		
	return(C)

def cofactor_matrix(M): 
	
	# M 4x4 matrix
	# return cofactor matrix of M

	C = np.empty(shape=[4, 4])
	
	for i in range(4):
		for j in range(4):
  			C[i][j] = ((-1)**(i+j)*determinant3(cofactor(M,i,j)))
	return(C)		

def adj_matrix(M): 

	# M 4x4 matrix
	# return adjungate matrix of M

	return(np.transpose(cofactor_matrix(M)))

def divisibility(M,q): 

	# M 4x4 integer matrix, q integer
	# return 1 if every entry of M is divisible by q
	# return 0 otherwise

	for i in range(4):
		for j in range(4):
			if (A[i][j] % q != 0):
				return(0)
	return (1) 

def permute(M,pe): 

	# M 4x4 matrix, pe 4x1 array - permutation of S4
	# return the matrix made from M by permutation its columns with permutation pe

	M2 = np.empty(shape=[4, 4])
	for i in range(4):
		for j in range(4):
			M2[i][j] = M[i][S4[pe][j]-1]
	return(M2)

def comparison(M,listL):

	# M 4x4 matrix, listL list of 4x4 matrices
	# return 0 if M is conatained in listL
	# return 1 otherwise

	l = len(listL)
	for k in range(l):
		if (subcomparison(M,listL,k) == 1):
			return(0)
	return(1)		

def subcomparison(M,listL,k):

	# M 4x4 matrix, listL list of 4x4 matrices, k integer (kth matrix in listL)
	# return 0 if kth matrix of listL is equal to M
	# return 1 otherwise

	for i in range(4):
		for j in range(4):
			if (M[i][j] != listL[k][i][j]):
				return(0)
	return(1)			

def text_output(listL,name):

	# listL list of 4x4 matrices, name string (name of the output textfile)
	# if listL is not empty, then print all matrices of listL to the textfile 'name'
	# otherwise print '0' to the textfile 'name'

	f = open(name,'w+')
	
	if (listL == []):
		f.write('0')
	else: 
		l = len(listL)
		for k in range(l):
			for i in range(4):
				for j in range(4):
					f.write('%d' % listL[k][i][j])
					f.write(" ")
				f.write("\n")
			f.write("\n")
	f.close()

#------------------------------------------------------------------------------------------
# MAIN PROGRAM

if ((first_checking(P) == 1) and (first_checking(Q) == 1)):
	list_of_solutions = []
	det_P = determinant4(P)
	adj_P = adj_matrix(P)
	det_Q = determinant4(Q)

	if ((det_P == det_Q) or (det_P == -det_Q)):	#otherwise P and Q cannot be lattice equivalent
		for perm in range(24):
			A = np.dot(permute(Q,perm),adj_P)
			if (divisibility(A,det_P)):
				A = A/det_P
				if (comparison(A,list_of_solutions) == 1):
					list_of_solutions.append(A)

	text_output(list_of_solutions,solution_file)
	print('The solution was printed to ',solution_file,'.')

else: print('Wrong input.')	
