#reset()

tau = (1 + sqrt(5))/2 #for working in SR
#tau = SR.symbol("tau", domain='positive') #for working in SR

def GenA(m): #matrix A con m columnas 
    mat = identity_matrix(m) +  block_matrix([[ 0, - identity_matrix(m-1) ], [ zero_matrix(1,1), 0 ]])
    #add first and last column
    mat = block_matrix(1, 3, [zero_matrix(m,1), mat, zero_matrix(m,1)], subdivide=False)
    mat[m-1,m+1]  = -1
    return mat

def MatGener(input,n):
    if (input in ['A','a'] and n>0):
        mat = GenA(n)
    elif (input in ['D','d'] and n>0):
        mat = block_matrix(2, 1, [GenA(n-1) ,zero_matrix(1,n+1)], subdivide=False)
        mat[n-1,n-1] = 1; mat[n-1,n] = 1;
    elif (input in ['E','e'] and n==8):
        mat = block_matrix(3, 1, [GenA(6) ,zero_matrix(1,8), ones_matrix(1,8)], subdivide=False)
        mat = block_matrix(1, 2, [ mat, zero_matrix(8,1)], subdivide=False)
        mat[6,6] = 1; mat[6,7] = 1; mat[6,0] = 0; mat[7,0] = 0; mat[7,8]= 1
    elif (input in ['F','f'] and n==4):
        mat = block_matrix(2, 1, [GenA(3) ,ones_matrix(1,5)], subdivide=False)
        mat[2,4] = 0; mat[3,0] = 0
    elif (input in ['H','h'] and n==3):
        mat = matrix([[0, 2, 0, 0], [0, -1-tau, 1, -tau], [0, 0, 0, 2]])
    elif (input in ['H','h'] and n==4):
        mat = block_matrix(2, 1, [ matrix([0, -1-tau, 1, 1, 1]), - GenA(3)], subdivide=False)
    else: 
        raise Exception('Error: not possible input')  
    return mat

def ListOfReflectors(mat):
    m = mat.ncols()
    listRef = []
    for row in mat:
        u = row/row.norm()
        Q = identity_matrix(m) - 2*(u.column()*matrix(u))
        listRef.append(Q)
    return listRef

def binary_insert(seq, t): 
    min = 0
    max = len(seq) - 1
    ind = -2
    while ind == -2:
        if max < min:
            ind = min
            Nlist = seq[0:ind] + [t] + seq[ind:] 
            return(Nlist, True)
        m = (min + max) // 2
        if seq[m] < t:    
            min = m + 1
        elif t < seq[m]:   
            max = m - 1
        else:
            ind = -1    
    return (seq, False)

def Orbit(X,n,v=[]):
    M = MatGener(X,n)
    ListRefl = ListOfReflectors(M)
    
    if v==[]:
        if X in ['a','A','e','E','D','d']:
            v =[i for i in reversed(range(M.ncols()))]
            v = matrix(v)
            v[0,0] = 1
        elif X in ['f','F']:
            v = matrix([ 1, 3, 2, 1, 8])
        else:
            if n==3:
                v = matrix(SR,[ 1, 1, 5, 1]) 
            else:    
                v = matrix(SR,[ 1, 0, 1, 2, 3]) 

               
    ##if X in ['a','A','e','E','D','d','h','H','f']:    # uncomment for testing what happens with f4
    if X in ['a','A','e','E','D','d','h','H']:          # comment for testing what happens with f4
        SetOfPoints = set([])
        count = 0
        CurrentList = [v]
        newPnt = True
        while (newPnt):
            newPnt = False
            NextSet = set([])
            for point in CurrentList:
                for matr in ListRefl:
                	#reflPnt = expand(matr*point.transpose()).simplify()  #use this instead of the following
                              # line when working with tau as symbolic (H3 and H4)
                              # commented because this code does not work with H3 nor H4
                	#reflPnt = matrix(SR, [ [reflPnt[i,0].simplify_full()] for i in range(reflPnt.nrows())]) 
                              # uncomment the previous line  when working with tau as symbolic  
                              # commented because this code does not work with H3 nor H4         
                    reflPnt = reflPnt.transpose()
                    reflPnt.set_immutable()
    
                    if reflPnt not in SetOfPoints:
                        ##if reflPnt == matrix([1, -4, -5, -6, 1]) : # uncomment for testing what happens with f4
                            ##print "Reflected point: %s"%(reflPnt)  # uncomment for testing what happens with f4
                            ##print "Is it already in the set of points? %s"%(reflPnt in SetOfPoints)#  ""
                        SetOfPoints.add(reflPnt)
                        NextSet.add(reflPnt)
                        count += 1
                        newPnt = True
                    reset('reflPnt')    
            NextSet = Set(NextSet) 
            CurrentList = NextSet.list()
            reset('NextSet')
           
        return count
              
    else: # Implemented for treating tau as float. 
          # It also works and gives the correct result for A_n, but it is slower.
          # The points in the orbit are stored ordered (lexigographycally) in lists of at most maxLenList elements.
          # Search inside these lists is done by binary search
          # ListMin stores the minimum element of those lists.
           
        maxLenList = 100
        ListPoints = [[],[v]] #list of sorted lists of points with less than maxLenList elements
        lenv=len(v.list())
        ListMin = [matrix(SR,zero_matrix(1,lenv)), v]
        count = 1
        CurrentList = [v]
        newPnt = True
        while (newPnt):
            newPnt = False
            NextList = []
            for point in CurrentList:
                for matr in ListRefl:
                    reflPnt = matr*point.transpose()
                    reflPnt = reflPnt.transpose()
                    
                    lenMin = len(ListMin)
                    for i in range(lenMin):
                        condi = False
                        if i==lenMin - 1: 
                            condi = True
                        else:
                            if reflPnt < ListMin[i+1]:
                                condi = True
                        if ((reflPnt >= ListMin[i]) and condi):
                            [ListPoints[i], ins] = binary_insert(ListPoints[i], reflPnt)
                            if ins:                               
                                NextList.append(reflPnt)
                                count += 1
                                newPnt = True
                                leng = len(ListPoints[i])
                                if leng >= maxLenList:
                                    indCut = (leng/2).ceil()
                                    ListPoints2 = ListPoints[0:i] + [ListPoints[i][0:indCut]] + [ListPoints[i][indCut:]] + ListPoints[i+1:] 
                                    ListPoints = ListPoints2                                    
                                    ListMin = ListMin[0:i+1] + [ ListPoints[i+1][0]] + ListMin[i+1:]
                            continue
            CurrentList = NextList                  
    return count

def check_orbit(filename):
    tests=[]
    with file(filename,'r') as f:
         cases = f.readlines()
         for line in cases:
             words = line.split()
             lenw = len(words)
             if lenw == 2:
                 currentTest = [words[0], int(words[1])] #it will give an error if the second word does not represent an integer
             else:
                 ini_vect = matrix([int(words[i]) for i in [2..(lenw - 1)]])
                 currentTest = [words[0], int(words[1]), ini_vect]
             tests.append(currentTest)
    for Test in tests:
        if len(Test) == 2:
            print "Test: %s%s"%(Test[0],Test[1])
            t = cputime()
            orb = Orbit(Test[0],Test[1])
        else:
            print "Test: %s%s; vector inicial = %s"%(Test[0],Test[1], Test[2])
            t = cputime()
            orb = Orbit(Test[0],Test[1],Test[2])
        tt = cputime(t)
        print "|G|: %s"%(orb)
        print "Time: %s"%(tt)
        print('')
    return None

check_orbit('OrbitTests.txt')

## uncomment the following lines for testing what happens with f4:
##orbF4_using_set = Orbit('f',4)
##print "Orbit obtained when using sets: %s"%(orbF4_using_set)
##orbF4_implemented_structure = Orbit('F',4)
##print "Orbit obtained when using the implemented structure: %s"%(orbF4_implemented_structure)

