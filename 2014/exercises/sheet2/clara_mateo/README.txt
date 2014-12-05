The functions reads to matrices from two given files (for instance, matrix_P1 and matrix_P2) and tries to find the matrix

X = (matrix_P1)^{-1} matrix_P2

Since we want det(X)=\pm 1 first it checks if abs(det(matrix_1))==abs(det(matrix_2)).
If it is, it computes the adjugate matrix of matrix_P1 = adj(matrix_P1) and check if the entries of

adj(matrix_P1) matrix_P2

are divisible by det(matrix_P1). In this case, we store the matrix X = 1/det adj(matrix_P1) matrix_P2.

Since the label of the points is important, we do the same for all the columns-permutations of matrix_P1.

Finally, we end up with a file called result_P1_P2.txt where all the transformations are saved in.