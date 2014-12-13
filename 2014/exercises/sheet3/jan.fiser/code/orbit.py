import numpy as np
import math
import time

def reflection(vector,normal_vector):
	''' Reflection of "vector" through a hyperplane defined by its "normal_vector". '''

	v = vector - np.multiply(2*(np.dot(normal_vector,vector)/np.dot(normal_vector,normal_vector)),normal_vector)
	return(v)

def generating_vectors(X,dim): 						#X_dim in Dn,An,E8,F4,H3,H4
	''' Returns list of generators of X_dim. '''

	generating_vectors = list()
	tau = 0.5 + math.sqrt(5)/2

	if (X == 'A'):
		for i in range(dim):
			a_0 = np.zeros(1+i, dtype=np.float)
			a_1 = np.array([1., -1.], dtype=np.float)
			a_2 = np.zeros(dim-i-1, dtype=np.float)
			a = np.concatenate((a_0,a_1,a_2), axis=0)
	
			generating_vectors.append(a)	
		return(generating_vectors)	

	if (X == 'D'):
		a_0 = np.zeros(dim-1, dtype=np.float)
		a_1 = np.array([1., 1.], dtype=np.float)
		a = np.concatenate((a_0,a_1), axis=0)

		generating_vectors.append(a)

		for i in range(dim-1):
			a_0 = np.zeros(1+i, dtype=np.float)
			a_1 = np.array([1., -1.], dtype=np.float)
			a_2 = np.zeros(dim-i-2, dtype=np.float)
			a = np.concatenate((a_0,a_1,a_2), axis=0)
	
			generating_vectors.append(a)

		return(generating_vectors)	

	if (X == 'E'):
		for i in range(6):
			a_0 = np.zeros(1+i, dtype=np.float)
			a_1 = np.array([1., -1.], dtype=np.float)
			a_2 = np.zeros(6-i, dtype=np.float)
			a = np.concatenate((a_0,a_1,a_2), axis=0)

			generating_vectors.append(a)	

		generating_vectors.append([0., 0., 0., 0., 0., 0., 1., 1., 0.])
		generating_vectors.append([0., 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
		
		return(generating_vectors)

	if (X == 'F'):
		generating_vectors.append([0.,  1.,  -1.,   0.,   0. ])
		generating_vectors.append([0.,  0.,   1.,  -1.,   0. ])
		generating_vectors.append([0.,  0.,   0.,   1.,   0. ])
		generating_vectors.append([0., -0.5, -0.5, -0.5, -0.5])
		
		return(generating_vectors)

	if (X == 'H'):
		if (dim == 3):			
			generating_vectors.append([0.,   2.,    0.,  0.])
			generating_vectors.append([0., -tau, 1/tau, -1.])
			generating_vectors.append([0.,   0.,    0.,  2.])

			return(generating_vectors)
		if (dim == 4):
			generating_vectors.append([0., -tau-1,  1.,  1., 1.])
			generating_vectors.append([0.,    -1.,  1.,  0., 0.])
			generating_vectors.append([0.,     0., -1.,  1., 0.])
			generating_vectors.append([0.,     0.,  0., -1., 1.])

			return(generating_vectors)

def vector_not_in_hyperplanes(X,dim):
	''' Returns a vector which doesn't lie in any of the reflecting hyperplanes of X_dim '''

	if (X == 'A'):
		vector_not_in_hyperplanes = np.zeros(dim+2, dtype=np.float)
	else: vector_not_in_hyperplanes = np.zeros(dim+1, dtype=np.float)

	for i in range(len(vector_not_in_hyperplanes)):
		vector_not_in_hyperplanes[i]=1.0 + i

	if ((X == 'F') or (X == 'E')): vector_not_in_hyperplanes[4] += 0.1

	return(vector_not_in_hyperplanes)

def array_to_string(array):
	''' Make a string of an array with float values. 
		Round the float values to given number of decimal places. '''

	array_to_string = ''
	for i in range(len(array)):
		array_to_string += str('%.10f' %array[i]) 		#here could be changed the number of decimal places
	return(array_to_string)	


# Reading input

input_string    = input('Xn = ')
group_letter 	= str(input_string[0])
group_dim 		= int(input_string[1])
input_vector 	= input('Optional input vector: ')

starting_time = time.time()

# Generators

generators = generating_vectors(group_letter,group_dim)

# Input vector

if (input_vector == ''):
	vector = vector_not_in_hyperplanes(group_letter, group_dim)
else: vector = np.array(input_vector.split(),dtype=np.float)

# Computing of the orbit

orbit = set()

orbit.add(array_to_string(vector))

stack = list()
stack.append(vector)

counter = 0

while(len(stack) > 0):
	v = stack.pop()
	for i in range(group_dim):
		reflected_v = reflection(v,generators[i])
		if (array_to_string(reflected_v) not in orbit):
			orbit.add(array_to_string(reflected_v))
			#counter +=1
			#print(counter)
			stack.append(reflected_v)		

# Output

print("The cardinility of the orbit is %s."%len(orbit))
end_time = time.time() - starting_time
print('Finished in %.2f seconds.'%end_time)