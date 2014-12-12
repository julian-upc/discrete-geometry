# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 21:53:32 2014

@author: Cassandra de Cock
"""

from matrix_generator import matrix_A, matrix_D, matrix_H, matrix_F
from vector_generator import vector_generator, space_generator, ref_tuple
import time


def orbit (X, n, v =[]):
    start_time= time.time()
    
    if X == 'A':
        M = matrix_A(n)
        n +=1
    elif X == 'D':
        M = matrix_D(n)
    elif X =='H':
        M = matrix_H(n)
    elif X == 'F' :
        M = matrix_F(n)
    M = list(M)

    if v:
        V = v
    else: 
        L = space_generator(M,n)
        M = L
        V = vector_generator(M,n)

    for i in range(len(V)):
        V[i] = round(V[i],6)

    ref_list = {tuple(V)}
    ref_list_new = {tuple(V)}
    while len(ref_list_new)!=0 : 
        ref_list_added = ref_list_new
        ref_list_new = set()
        for v in ref_list_added:

            for a in M:                
                ref_v = ref_tuple(v,a)
                
                ref_list_new = ref_list_new.union({ref_v}-ref_list)
                ref_list = ref_list.union({ref_v})
                
    print(len(ref_list))
    print(time.time()-start_time)
