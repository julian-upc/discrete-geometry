compile with the order: g++ reflections2.cc -o reflections2.out

execute: ./reflections2.out
         to see time performance: time ./reflections2.out (look at user time).

Input:  one uppercase letter A,D,E,F,H.
        integer number as dimension n>=2 for A and D, n=8 for E, n=4 for F, n=3,4 for H.
	either 0 or one vector of dimension n+1 (n+2 for A) that does not start with 0.

If the input is not valid an error message shoud appear.

Outputs: The seize of the orbit the vector given through the reflection group selected or the seize of the orbit of one vector that does not lie in any of the reflecting hyperplanes of the given reflection group.

Tests: It works and it is fast for A, D <8, F 4 and H 3. It does not work due to precision issues for H 4, and it has not been tested due to velocity issues for E 8, and A, D >8.
