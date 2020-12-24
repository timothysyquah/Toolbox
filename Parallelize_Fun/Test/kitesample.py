#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 16:19:47 2020

@author: tquah
"""


import multiprocessing
import time as time
import math as m
import numpy as np
import matplotlib.pyplot as plt

n = np.logspace(1,7.5,14,dtype = int)
repeat = 4
pool = multiprocessing.Pool(4)
result = np.zeros((len(n),2))
# result = np.zeros((len(n)))



def fun(i):
    return i*i*i/i/i



def nonparallel(n):
    a1 = []
    #nonparallel method
    start = time.time()
    for i in range(1,n,1):
        a1.append(fun(i))
    end = time.time()
    return end-start
    
def parallel(n):
    start = time.time()
    a2= pool.map(fun, range(1,n))
    end = time.time()
    return end-start

method = [nonparallel,parallel]

for i in range(0,len(n),1):
    
    
    for j in range(0,len(method),1):
        temp  = []
        for k in range(0,repeat):
            t = method[j](n[i])
            if k!=0:
                temp.append(t)
        result[i,j] = np.mean(temp)
plt.close('all')
plt.figure()
plt.scatter(np.log(n),np.log(result[:,0]),label = 'Serial')
plt.scatter(np.log(n),np.log(result[:,1]),label = 'Parallel 4 Cores')
plt.xlabel('$\log(n)$')
plt.ylabel('$\log(t)$')
plt.legend()




# a1sort = sorted(a1)
# a2sort = sorted(a2)

# if a1sort!=a2sort:
#     print('Output Not Same!')