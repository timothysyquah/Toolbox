#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 17:17:15 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt 

def fun(x,a):
    c = x+a**2
    print(c)
    d = c+1
    return d*x

a = np.linspace(0,10,10)
x = np.linspace(0,10,100)
plt.figure()

for i in range(0,len(a),1):
    y = fun(x,a[i])

    plt.plot(x,y)