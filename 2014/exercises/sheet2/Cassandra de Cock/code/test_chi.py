# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 20:34:51 2014

@author: Cassandra de Cock de Rameyen
"""
import numpy as np
from chi import chi

def test_chi(matrix_list):
    
    list_matrix = np.loadtxt(matrix_list, dtype='str')
    size = list_matrix.shape
    for i in range(size[1]):
        chi(list_matrix[i,0],list_matrix[i,1],list_matrix[i,2])
    
    