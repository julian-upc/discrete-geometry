# -*- coding: utf-8 -*-
"""
Created on Wed Dec  3 22:07:02 2014

@author: Cassandra de Cock
"""
import numpy as np

def dot_tuple(v,l):
    temp = 0
    for i in range(len(v)):
        temp += v[i]*l[i]
        
    return temp

def ref_tuple(v,a):
    temp1 = dot_tuple(v,a)
    temp2 = dot_tuple(a,a)
    
    ref_v = []  
    for i in range(len(v)):
        ref_v.append(round(v[i]-2*temp1/temp2*a[i],6))

    return(tuple(ref_v))


def space_generator(M,n):
    
    gen_list = set()
    for i in range(len(M)):
        gen_list = gen_list.union({tuple(M[i])})
        
    for i in range(len(M)):
        gen_list_new = {tuple(M[i])}

        while gen_list_new : 
            gen_list_added = gen_list_new
            gen_list_new = set()
            
            for v in gen_list_added:
                for i in range(len(M)):
                    a = M[i]
                    
                    ref_v = ref_tuple(v,a)
                    
                    gen_list_new = gen_list_new.union({ref_v}-gen_list)
                    gen_list = gen_list.union({ref_v})

    return(gen_list)
 
def vector_generator(L, n):

    v = np.random.randint(2*n, size = n+1)
    v[0]=1
    
    not_on_one = True
    
    for l in L:
        if dot_tuple(v,l) == 0:
            not_on_one = False
        
    while not_on_one == False : 
        v = np.random.randint(2*n, size = n+1)
        v[0]=1
        not_on_one = True
    
        for l in L:
            if dot_tuple(v,l) == 0:
                not_on_one = False
        
    return v
