#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:34:30 2020

@author: tquah
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
path3_dob = '/home/tquah/Projects/PhaseDiagram/Dobreninpts.dat'
dobarray = np.loadtxt(path3_dob)

y = np.sqrt(np.power(dobarray[:,1],4))
x = 1-dobarray[:,0]

plt.close('all')
plt.figure()
plt.plot(x[:35],y[:35],'--r')
plt.plot(x[36:77],y[36:77],'--r')
plt.plot(x[78:134],y[78:134],'--r')
plt.plot(x[135:],y[135:],'--r')

# plt.text(0.025,2,'S',color = 'r',size = textsize)
# plt.text(0.18,2.5,'C',color = 'r',size = textsize)
# plt.text(0.4,3.5,'L',color = 'r',size = textsize)


pathleq = '/home/tquah/Projects/PhaseDiagram/jleq.csv'
df = pd.read_csv(pathleq)
items = list(df)

for i in range(0,int(len(items)/2)):
    xraw = np.array(df[items[2*i]])
    x = xraw[~np.isnan(xraw)]
    yraw = np.array(df[items[2*i+1]])

    y = yraw[~np.isnan(yraw)]
    plt.plot(x,y,'k')
plt.ylim(1,4.5)
plt.xlim(0,0.6)

# leqarray = np.loadtxt(pathleq,delimiter = ',')