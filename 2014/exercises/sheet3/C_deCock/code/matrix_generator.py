# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 21:34:36 2014

@author:Cassandra de Cock
"""
import numpy as np
import math

def matrix_A(n):
    A = np.zeros((n,n+2))
    for i in range(n):
        A[i,i+1] = 1
        A[i,i+2] = -1
    
    return A

def matrix_D(n):
    D = np.zeros((n,n+1))
    
    for i in range(n-1):
        D[i,i+1] = 1
        D[i,i+2] = -1
    
    D[n-1,n-1]=  1
    D[n-1,n] = 1
    return D
    
def matrix_H(n):
    tau = 0.5 + 0.5*round(math.sqrt(5),6)
    if n==3:
        H = np.array([[0,2,0,0], [0, -tau**2, 1, -tau], [0,0,0,2]])
        return H
    elif n==4 :
        H = np.array([[0, -tau**2, 1, 1, 1], [0, -1, 1, 0,0], [0,0,-1,1,0],[0,0,0,-1,1]])
        return H
    else :
        print('Matrix not known')

    
def matrix_F(n):
    if n ==4 :
        F = np.array([[0, 1, -1, 0, 0], [0, 0, 1, -1, 0], [0, 0, 0, 1, 0], [0, -0.5, -0.5, -0.5, -0.5]])
        return F
    else:
         print('Matrix not known')
         
    