#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  2 14:36:20 2020

@author: tquah
"""

import numpy as np
from sklearn.linear_model import LinearRegression
from patsy import cr
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def powerlaw(x,gamma):
    return np.power(x,gamma)

def powerlaw_random(x,gamma):
    return np.power(x ,gamma)+0.0001*np.random.normal(0, 1,len(x))

plt.close('all')
x = np.linspace(0,0.5,100)
y = powerlaw(x, 2/3)
y_noise = powerlaw_random(x, 2/3)

param,pcov = curve_fit(powerlaw,x,y_noise)



plt.figure()
plt.plot(x,y,'k')
plt.scatter(x,y_noise,color = 'r',alpha = 0.5)
    
n_obs = len(x)

df = 10
# Generate spline basis with 10 degrees of freedom
x_basis = cr(x, df=df, constraints='center')
# Fit model to the data
model = LinearRegression().fit(x_basis, y_noise)
# Get estimates
y_hat = model.predict(x_basis)

# plt.scatter(x, y)
plt.plot(x, y_hat, 'r')
plt.title(f'Natural cubic spline with {df} degrees of freedom')
