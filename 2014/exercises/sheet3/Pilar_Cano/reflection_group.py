#!/usr/bin/python

from numpy import *
from sets import *
from numpy.random import * 

#This function checks if the vector is in the planes given with the normal vectors in M, where M is a set and v the vector.
def Coplanar(M,v):
   y=1
   for x in M:
     d = float(dot(array(x),array(v)))
     if d == 0:
       y = 0
   return (y)
   
   
#It suppose to calculate the orbit of a vector given X=A,D,E,F or H.
def Orbit(X,n,v):
   tau=float(1)/2+(float(1)/2)*((5)**(.5))
   k=n+1
   
   #Is generating An
   if X == 'A':
      k=n+2
      M=zeros((n,n+2))
      for i in range (0,n):
        for j in range(1,n+1):
          if j==i+1:
            M[i][j]=1
            M[i][j+1]=-1
    
   #Is generating Dn        
   if X == 'D':
      M=zeros((n,n+1))
      for i in range (0,n-1):
        for j in range(1,n):
          if j==i+1:
            M[i][j]=1
            M[i][j+1]=-1
        M[n-1][n-1]=1
        M[n-1][n]=1
    
   #Is generating F4 
   if X == 'F':
      M = [[0,1,-1,0,0],
           [0,0,1,-1,0],
           [0,0,0,1,0],
           [0,-0.5,-0.5,-0.5,-0.5]] 
   
   #Is generating Hn
   if X == 'H':
      if n == 3:
         M = [[0,2,0,0],
              [0,-tau,float(1)/tau,-1],
              [0,0,0,2]]
      if n == 4:
         M = [[0, -tau , float(1)/tau, float(1)/tau, float(1)/tau],
              [0, -1, 1, 0, 0],
              [0, 0, -1, 1, 0],
              [0, 0, 0, -1, 1]]
              
   #The generate matrices for each En are the ones given in wikipedia
   if X == 'E':
     if n == 6:
        M = [[0,1,-1,0,0,0,0],
             [0,0,1,-1,0,0,0],
             [0,0,0,1,-1,0,0],
             [0,0,0,0,1,-1,0],
             [0,0,0,0,1,1,0],
             [0,-0.5,-0.5,-0.5,-0.5,-0.5,float(3**(.5))/2]]
     if n == 7:
        M = [[0,1,-1,0,0,0,0],
             [0,0,1,-1,0,0,0],
             [0,0,0,1,-1,0,0],
             [0,0,0,0,1,-1,0],
             [0,0,0,0,1,1,0],
             [0,-0.5,-0.5,-0.5,-0.5,-0.5,float(1)/((2)**(.5))]]
     if n == 8:
        M = [[0,1,-1,0,0,0,0],
             [0,0,1,-1,0,0,0],
             [0,0,0,1,-1,0,0],
             [0,0,0,0,1,-1,0],
             [0,0,0,0,1,1,0],
             [0,-0.5,-0.5,-0.5,-0.5,-0.5,-0.5]]
     
   
   #Is finding a vector that it isn't in the reflecting hyperplanes        
   if v==[]: 
     ref_hp=set([])
     for i in range(n):
       ref_hp = ref_hp|set([tuple(M[i])])
     S = ref_hp
     while (S != set()):
       for i in range (0,n):
         for x in S:
           g = array(x)+((-2)*((float((dot(array(x),array(M[i]))))/(dot(M[i],M[i])))*array(M[i])))
           g = around(g, decimals=6)
           S= S|set([tuple(g)])
       S = S - ref_hp
       ref_hp = ref_hp|S
     y = 0
     while(y == 0):
       v = randint(2*n, size=k)
       y = Coplanar(ref_hp,v)
    
   #Is making all the reflections of v  
   t=tuple(v)
   ref=set([t])
   ref_new=set([t])
   while (ref_new != set()):
     for i in range (0,n):
       for x in ref_new:
         r = array(x)+((-2)*((float((dot(array(x),array(M[i]))))/(dot(M[i],M[i])))*array(M[i])))
         r = around(r, decimals=7)
         l = tuple(r)
         ref_new=ref_new|set([l])
     ref_new=ref_new - ref
     ref = ref|ref_new
     o = len(ref)
     print(o)
     
  

       