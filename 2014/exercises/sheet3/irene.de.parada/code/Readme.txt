Function check_orbit is given the name of a file containing in each row one test. 
 The first character in each test must be the letter of the group (it can be lowercase
 or uppercase) and the second one the number corresponding to the dimension, n. 
 Optionally we can write after in the same line the numbers corresponding to the 
 homogeneous coordinates of the starting vector (otherwise one producing a maximum
 cardinality orbit will be taken). 
 
Function Orbit(X,n,[v_inicial]), where X is the letter lowercase or uppercase, is then
 called. If Xn corresponds neither to An, Dn, E8, F4, H3 nor H4, an exception is raised.
 However, even if an exception is not raised, this code cannot deal with H3 nor H4.
 
 In each case the reflection matrices corresponding to the generating reflecting 
 hyperplanes are computed and stored. In each step, all the new points obtained in the
 previous step (we begin with v_inicial) are reflected in all that hyperplanes. 
 
 All different points obtained are stored in different structures depending on the case.
 If Xn corresponds to An, Dn or E8, a set is used for this purpose. Otherwise, a ordered
 list of ordered lists with a maximum number of elements is used. Binary search is used 
 inside. 
 
 With F4 the second structure is used because the first doesn't work. For seeing what
 happens when using sets there are some commented lines in the code with instructions
 that can be uncommented. There is also a commented test for this case at the end of the
 code.
 
 The code is designed for working with the tau in generating matrices of H3 and H4 as 
 symbolic. However, it doesn't work in these cases, so the lines corresponding to 
 simplifying the expressions with symbolic elements are commented. 


** Instructions
Run the .script from the command line using 'sage Orbit.sage'. 

Upon start, the program will read the required tests from a file called 'OrbitTests.txt' 
which must be located on the same folder.