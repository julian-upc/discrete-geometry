In order to execute tetrahedron.py take into account that:
instructions:
-open python
-import tetrahedron
-from tetrahedron import main
-execute main with the following command:
	-main(‘yourfileP.txt’,’yourfileQ.txt’,’outputfilename.txt’)
where fileP.txt and fileQ.txt have to be files with the input 4x4 matrices made up of integers which contain 4 points in general position. If the inputs are not ok an error message will be printed either in outputfile.txt or as an error of a function in the terminal.
-If the inputs are given correctly as in matrix_P.txt or matrix_Q.txt the output file will be a txt with all the possible A for the sets of points given. If no matrix A exist a 4x4 matrix of 0’s will be given as output instead.


In order to execute testchi.py take into account that:

-testchi does not check the inputs, as it is not intended for a user but for a personal use. Three files with matrices must be given, two of them with 4 points P and Q and another matrix A with the Z affine map. 

 There is a very basic main that given the inputs says Yes! if A transforms P to Q or No!if not and a function testchi that can be called from another function or main etc that given the 3 files returns true or false.

The other files extension .txt are either matrix inputs or outputs that have been used to test the code.