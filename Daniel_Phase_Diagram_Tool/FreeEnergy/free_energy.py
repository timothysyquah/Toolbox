#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 22:45:21 2020

@author: tquah
"""
#https://stackoverflow.com/questions/51321100/python-natural-smoothing-splines
import numpy as np
from sklearn.linear_model import LinearRegression,RANSACRegressor
from patsy import cr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
plt.close('all')
import itertools

def stripwhitespace(lst):
    lst = [x.strip('') for x in lst]
    indices = [i for i, x in enumerate(lst) if x == '']
    for index in range(len(indices)-1,-1,-1):
        lst.pop(indices[index])
    return lst
def powerlaw(x,a,b,c,d):
    return a*np.power(x,b)+c*x+d



plt.close('all')
path = 'free_energy.dat'
op = open(path,'r')
rawdata = op.read().split('#')
op.close()
# marker = itertools.cycle(('+', 'o', 's','v','^','<','>','x','p','h','H','*'))
c  = ['b','r','g','k','m','c']


data_dictionary = dict()
for i in range(2,len(rawdata),1):
    raw = rawdata[i].split('\n')
    raw = stripwhitespace(raw)
    array = np.zeros((len(raw)-1,2))
    for j in range(1,len(raw),1):
        # print(raw[j].split(' '))
        temparray =  np.array(raw[j].split(' '),dtype=float)
        
        array[j-1,0]=  temparray[0]
        array[j-1,1]=  temparray[1]
    indices = np.argsort(array[:,0],0)
    array_sorted = array[indices,:]        
    data_dictionary[raw[0]] = array_sorted
    
    
    
    
    
    
    
    
    
    
    
    
    
# plt.scatter(array_sorted[:,0],array_sorted[:,1])
# param,pcov = least_squares(powerlaw,array_sorted[:,0],array_sorted[:,1])
# newy = powerlaw(array_sorted[:,0],param[0],param[1],param[2],param[3])

header = list(data_dictionary)
DIS_index = header.index('DIS')
DIS_array = data_dictionary[header[DIS_index]]
plt.figure()


phaselist = []

for i in range(0,len(header),1):
    if header[i]!='DIS':
        print(header[i])
        xx = data_dictionary[header[i]][:,0]
        # print(x)
        index = np.intersect1d(DIS_array[:,0],xx,return_indices=True) 
        
        
        yy = data_dictionary[header[i]][:,1]-DIS_array[index[1],1]

        loc = np.where(np.abs(yy)>1e-10)[0]
        x = xx[loc]
        y = data_dictionary[header[i]][loc,1]
        df = len(x)-1
        # Generate spline basis with 10 degrees of freedom
        x_basis = cr(x, df=df, constraints='center')
        # Fit model to the data
        model = RANSACRegressor().fit(x_basis, y)
        # Get estimates
        y_hat = model.predict(x_basis)
        
        # plt.scatter(x, y)
        # plt.plot(array_sorted[:,0], y_hat, c[i])
        xplot = np.linspace(np.min(x),np.max(x),1000)
        x_basis_large = cr(xplot, df=df, constraints='center')
        y_hat = model.predict(x_basis_large)
        
        cs = CubicSpline(x,y)
        yplot = cs(xplot)
    # # Compare estimated coefficients

        # plt.plot(xplot,y_hat,color = c[i])
        plt.plot(xplot,yplot,'--',color = c[i])

        plt.scatter(x,y,marker='>',color = c[i],label = header[i])
        # plt.scatter(array_sorted[:,1],newy)
        phaselist.append(cs)
    # # else:
plt.legend()
    
    
    
    #     df = 4
    #     # Generate spline basis with 10 degrees of freedom
    #     x_basis = cr(array_sorted[:,0], df=df, constraints='center')
    #     # Fit model to the data
    #     model = LinearRegression().fit(x_basis, array_sorted[:,1])
    #     # Get estimates
    #     y_hat = model.predict(x_basis)
        
    #     # plt.scatter(x, y)
    #     # plt.plot(array_sorted[:,0], y_hat, c[i])
    #     xplot = np.linspace(np.min(array_sorted[:,0]),np.max(array_sorted[:,0]),100)
    #     x_basis_large = cr(xplot, df=df, constraints='center')
    #     y_hat = model.predict(x_basis_large)
        
        
    #     plt.plot(xplot,y_hat,color = c[i])
    #     plt.scatter(array_sorted[:,0],array_sorted[:,1],marker='>',color = c[i])
    #     # plt.scatter(array_sorted[:,1],newy)
    
    
    