#!/usr/bin/python

from numpy import *
from itertools import permutations


#Here I am permutating the columns of my matrix.
def Perm(A):
  V=[0, 1, 2, 3]
  T=[]
  P=list(permutations(V))
  for p in P:
    A1 = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]
    s = list(p)
    for i in range(0, 4):
      for j in range(0, 4):
        A1[i][j] = A[i][s[j]]
    T.append(A1)
  return T 
  
#This function calculate the determinant of a cuadratic matrix with size n
def Det(A, n):
  d = 0
  if n == 2:
    d = ((A[0][0])*(A[1][1]))-((A[0][1])*(A[1][0]))
    return d
  #Calculating the cofactors
  for i in range(0,n):
    M = delete(delete(array(A),i,0),0,1)
    q = Det(M,n-1)
    t = (((-1)**(i))*(A[i][0]))*(q)
    d = d + t
  return d

#This function heck if two cuadratic matrices with same size are the same 
def ComMa(A,T,n):
  for i in range(0,n):
    for j in range(0,n):
      if A[i][j] != T[i][j]:
        return 0
  return 1

#This function check if A is in M  
def ElemSet(A,M):
  m = len(M)
  n = len(A)
  if m == 0:
    return 1
  for i in range(0,m):
    y = ComMa(A,M[i],n)
    if y == 1:
      return 0
  return 1
 
#This funtion check if the entries of a square matrix are divisble by d  
def Div(A,d):
  m = len(A)
  for i in range(m):
    for j in range(m):
      w = abs(A[i][j]%d)
      if w != 0:
        return 0 
  return 1 
 
#This check if the entries of a square matrix are integers  
def Integ(A):
  m = len(A)
  for i in range(m):
    for j in range(m):
      if A[i][j] != int(A[i][j]):
        return 0
  return 1
            

#This check if two tetrahedra are Lattice Equivalent
def LatEq(T, Q, R):
   f = open(R, 'w')  
   Pt = transpose(T)
   Qt = transpose(Q)
   
   Z = [[0,0,0,0],
        [0,0,0,0],
        [0,0,0,0],
        [0,0,0,0]]
   e = Integ(T)
   e2 = Integ(Q)
   if e != 1:
     print('I only work with lattice tetrahedron.')
     f.write('I only work with lattice tetrahedron.')
     return False
   if e2 != 1:
     print('I only work with lattice tetrahedron.')
     f.write('I only work with lattice tetrahedron.')
     return False
   #Check if P and Q are tetrahedra and not polygons
   dp = Det(Pt, 4)
   dq = Det(Qt, 4)
   if dp == 0:
    print('I only work with lattice tetrahedron.')
    f.write('I only work with lattice tetrahedron.')
    return False
   if dq == 0:
    print('I only work with lattice tetrahedron.')
    f.write('I only work with lattice tetrahedron.')
    return False
    
   #If the determinant is different, then they cannot be equivalent
   if abs(dp) != abs(dq):
     f.write(str(matrix(Z)))
     print Z
     return False
   
   #Do the permutations with it respective A  
   L=Perm(Pt)
   A=[]
   for l in L: 
     Pad = [[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0]]
     dl = Det(l, 4)
     for i in range(0,4):
       for j in range(0,4):
         M = delete(delete(l,i,0),j,1)
         x = Det(M , 3)
         Pad[i][j] = (((-1)**(i+j))*(x))       
     Padt = transpose(Pad)  
     Al = dot(Qt,Padt)
     h = Div(Al,dl)
     if h == 1:
       Anew = Al/dl
       w = Det(Anew, 4) 
     #Check if the det of Al is +-1
       if abs(w) == 1:
     #check if Al is in A
         y = ElemSet(Anew, A)
         if y == 1: 
           print(Anew)   
           A.append(Anew)
           f.write('\n')
           f.write(str(Anew))
           f.write('\n')
           #printMatToFile(Anew, f)
     #printMatToFile(A,f)        
   m = len(A)
   #Check if P and Q are Lattice Equivalent
   if m == 0:
     print(Z)
     f.write(str(matrix(Z)))
   f.close()    
   #return A
   
  
  
def TestSuit(T,M,R):
  Q = loadtxt(M)
  P = loadtxt(T)
  Y = LatEq(P, Q, R)


  
  
  
print('Hello')  

TestSuit('Matrix_P1.txt','matrix_P3.txt','result_P1_P3.txt')  
  
  
  
   
   
   

P=[[1, 0, 0, 0],[1, 1, 0, 0], [1, 0, 1, 0], [1, 2, 3, 5]]   

Q=[[1, 0, 0, 0], [1, 0, 1, 0], [1, 1, 0, 0], [1, 2, 3, 5]]   
            
F =[[1, 0, 0, 0],[1, -1, 0, 0],[1, 0, 1, 0],[1, 2, 3, 5]]   

T =[[1, 0, 0, 0], [3, 0, 0, 0], [1, 0, 1, 0],[1, 2, 3, 5]] 

G = [[1,0,0,0], [1, 3, 0, 0], [1, 0, 1, 0], [1, 2, 3, 5]]

R =[[1, 0, 0, 0],[1, -1, 0, 0],[1, 0, 1, 0],[1, -2, 3, 5]]

