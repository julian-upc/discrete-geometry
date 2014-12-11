              orbit.py

          Jan Fiser, 11/12/2014

----------------------------------------

TASK

The program orbit.py computes the cardinality of the orbit of a vector v under X_n,
where X_n is one of the reflection groups:
A_n, D_n, E_8, F_4, H_3, H_4.

----------------------------------------

USER MANUAL

Run orbit.py.

On the output is displayed 'Xn = '.
Write there one of the following: 
An, Dn, E8, F4, H3 or H4 (where n is a natural number) 
and press Enter.

On the output is displayed 'Optional input vector: '.
Either press Enter to choose the predefined vector, 
whose orbit has cardinality |X_n|,
or write there an vector (of the right length), e.g.
5.123 6 1.0012 52. 8.1, 
and press Enter.

On the output is printed the cardinality of the orbit and the running time.

----------------------------------------

REMARKS

1/ This program cannot handle H_4 in most cases (also when the predefined vector is used).

2/ The vector which doesn't lie in any of the reflection hyperplanes of X_n
was chosen as [1 2 3 .. n+1] (in F4 and E8 as [1 2 3 4 5.1] and [1 .. 4 5.1 6 .. 9]).
These vectors satisfy the given condition (proven experimentally).

3/ I wanted to store all computed elements of the orbit in a set. 
But the data structure "set" cannot contain arrays (which corresponds to vectors). Therefore, there are stored strings representing the arrays in the set "orbit".
The string is made as follows:
	1. round all float values in the array to a given number of decimal places,
	2. consider all these decimal values as strings and concatenate them.

4/ There are no controls of the input in the code. The program expects a correct input
(e. g. as an input vector a vector of length n+1 (or n+2 for the group An) with int or float values). 
