#reset()

def stringRepresentsInteger(input):
    charNumbers = [str(i) for i in range (10)]
    if input[0]=='-':
        input=input[1:]
    for c in input:
        if not(c in charNumbers):
            return False
    return True

def columnIsZero(mat):
    lzero = vector([0] * mat.nrows()).list()
    for i in range( mat.ncols() ):
        if mat[:,i].list() == lzero:
            return True
    return False

def transfHomogeneous(matriz):
    mat = copy(matriz)
    for col in range( matriz.ncols() ):
        if ((matriz[0,col] != 0) and (matriz[0,col] != 1)):
            colu = vector(matriz[:,col].list())
            for row in range(matriz.nrows()):
                if ((colu[row] % matriz[0,col]) != 0):
                    raise Exception('Error reading matrix: it must define a lattice polytope')             
                mat[row,col] = round(colu[row]/matriz[0,col])          
    return mat

def readCheck(filename): #reads matrix checking if it is valid
    mat = []  
    with file(filename,'r') as f:
        for line in f:
            currentRow = []
            for element in line.split():
                if not stringRepresentsInteger(element):
                    raise Exception('Error reading matrix: elements must be integers')    
                currentRow.append(Integer(element))
            if len(currentRow) != 4 or len(mat) == 4:
                raise Exception('Error reading matrix: matrix must be 4x4') 
            mat.append(currentRow) 
        if len(mat) != 4:
             raise Exception('Error reading matrix: matrix must be 4x4') 
    matr = matrix(mat).transpose()         
    if columnIsZero(matr):
        raise Exception('Error reading matrix: column is zero vector')                        
    return matr

def matrixOfIntegers(mat):
    for row in mat:
        for element in row:
            if not stringRepresentsInteger(element.str()):
                return False
    return True

def permuMatrix(mat): #list with all possible matrices with the columns permuted
    matT = copy(mat).transpose() 
    listOfMatrices = []
    for perm in Permutations([0,1,2,3]):
        listI = [matT[perm[i]] for i in range(4)]
        listOfMatrices.append(matrix(listI).transpose())
    return listOfMatrices

def writeMatrixValues4(m, filename):
     for row in m.rows():
         #print("%d %d %d %d" % (row[0], row[1], row[2], row[3]))
         rowString = `row[0]`
         for i in range(1, len(row)):
            rowString += ' ' 
            rowString += `row[i]`
         rowString += '\n'
         filename.write(rowString)

def test_chi(filename):
    tests = [] 
    solutionList = [] 
    with file(filename,'r') as f:
        for line in f:
            currentTest = []
            for word in line.split():
                currentTest.append(word)
            if len(currentTest) == 3:
               tests.append(currentTest)
    counter = 0
    for currentTest in tests:
        resultFile = file(currentTest[2],'w')
        counter += 1
        #resultFile.write('\n \n \n' + str(counter) + ''.join('-' * 20 ) + '\n')
        P = readCheck(currentTest[0])
        Q = readCheck(currentTest[1])
        P = transfHomogeneous(P)
        Q = transfHomogeneous(Q)
        if (rank(P)<4 or rank(Q)<4):
            raise Exception('Error reading matrix: matrix must define a 3-dimensional tetrahedron but points lie in same plane') 
        possibleMatrices=set([])    
        for mP in permuMatrix(P):
            AA = Q*(mP.inverse()) 
            A = AA.delete_columns([0]).delete_rows([0])
            if (matrixOfIntegers(AA) and ([1,0,0,0] == AA[0,:].list()) and (abs(A.determinant())==1)):
                AA.set_immutable()
                possibleMatrices.add(AA)
                reset('AA')
        possibleMatrices = Set(possibleMatrices)        
        for mat in possibleMatrices.list():
             writeMatrixValues4(mat, resultFile)
             resultFile.write('\n')
        solutionList.append(possibleMatrices.list())
        reset('possibleMatrices')
        resultFile.close()
        #print('-----------------------')
        #print(possibleMatrices)
        #print('-----------------------')
    
    return solutionList

tests = test_chi("test_chi_input.txt")

#print(tests)
