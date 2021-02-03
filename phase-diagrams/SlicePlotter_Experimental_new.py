#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 29 13:49:23 2021

@author: tquah
"""

from matplotlib.lines import Line2D

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] =23
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 11
plt.close('all')

plt.close('all')
plt.figure(figsize=(7,5))
yspace = np.linspace(0,1.0,100)
data_array = np.loadtxt('asymslice.dat')
for i in range(0,len(data_array[:,0])):
    plt.plot(data_array[i,0]*np.ones_like(yspace),yspace,'--k',alpha = 0.3)
plt.yticks(color='w')
plt.xlabel('$f_A$')
overallphases = ['L','C','S','D']

path_exp = '/media/tquah/TOSHIBA EXT/Projects/PhaseDiagram/Experimental_Dataset_v2.csv'

df = pd.read_csv(path_exp)
header = list(df)

references = sorted(list(set(list(df['Reference']))))
update_references = []
for ref in references:
    newref = ref.replace('-',' ')
    update_references.append(newref)
               
fullreference = list(df['Reference'])
fullphases =  list(df['Phase'])
                  
                  
phases = sorted(list(set(list(df['Phase']))))
#symbol = ['h','+','+','h']
symbol = ['x','>','+','o']

color = ['k','y','g','r','b']

ab_ratio = np.sqrt(np.array(df['NscB']/df['NscA']))*np.array((df['bB']/df['bA']))

NscAverage =  np.array((df['NscB']+df['NscA'])/2)

Nsc_Backbone = NscAverage/np.array(df['Nbb(Total)'])

NaNloc = np.where(np.isnan(Nsc_Backbone)==True)[0]
Nsc_Backbone[NaNloc]=0


count1 = 0
x = np.arange(10)
condlist = [x<3, x>5]
choicelist = [x, x**2]
# print(np.select(condlist, choicelist))





for ref in references:
    count2 = 0 
    for phase in phases:
#        print(phase)
        refindex= [i for i, x in enumerate(fullreference) if x ==ref]
        phaseindex= [i for i, x in enumerate(fullphases) if x ==phase]
        height = overallphases.index(phase)
        loc = list(set(refindex) & set(phaseindex))
        if len(loc)>0:
            
            # print(Nsc_Backbone[loc])
            
            

            mloc = np.where(Nsc_Backbone[loc]<5.0)[0]
            floc = []
            for m in mloc:
                floc.append(np.where(Nsc_Backbone[loc][m]==Nsc_Backbone)[0][0])
            if len(floc)>0:
                
                yraw = ab_ratio[floc]
                # print(yraw)
                xraw = np.array(df['fA'][floc])
    
                yloc = np.where(yraw>1)[0]
                xfilter = xraw[yloc]
                yfilter = yraw[yloc]
                
#                print(ref)
#                print(floc)

                
                
                # print(yfilter)
                # nscnbbloc = np.where(Nsc_Backbone[loc]<1.0)[0]
                # print(nscnbbloc)
            
                # print(xfilter)
                name = phase+'-'+ref
                if len(xfilter)>0:
                    print((height/30+0.2))

                    plt.scatter(xfilter,(height/30+0.2)*np.ones_like(xfilter),c=color[count2],marker = symbol[count2],label = name,s = 100)
                
        count2+=1
    count1+=1

ypos = 6
# plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
# plt.savefig('/home/tquah/Figures/slice_otherdata.pdf',dpi = 300)

# plt.figure()

# for i in range(0,len(data_array[:,0])):
#     plt.plot(data_array[i,0]*np.ones_like(yspace),yspace,'k')
plt.yticks(color='w')
plt.xlabel('$f_A$')
# arraylist  = [np.array([1-0.384,1-0.457]),\
#                np.array([0.241]),\
#                    np.array([0.3,1-0.651,1-0.716,1-0.771]),\
#                       np.array([1-0.883])]
#plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

xmin = 0.08
xmax = 0.92
xarray = np.ones(10)
xarray[0] = xmin
sortdata = np.sort(data_array[:,0])
xarray[1:-1] = sortdata
xarray[-1] = xmax
textsize = 20
order = ['D','S','C','G','L','G','C','S','D']
avg = lambda x: (x[1]+x[0])/2

y= 0.15
for i in range(0,len(xarray)-1,1):
    xpos = avg(xarray[i:i+2])
    plt.text(xpos,y,order[i],color = 'k',size = textsize,ha = 'center',alpha = 0.7)

def plotbracket(x1,x2,y,h,label):
    plt.plot([x1, x1, x2, x2], [y-h, y, y, y-h],'-k', lw=1.5)
    plt.text((x1+x2)*.5, y, label, ha='center', va='bottom',size = textsize, color=col)

plt.xlim(xmin,xmax)
col = 'k'
h = 0.05
y = 0.3
x1, x2 = 0.105, 0.13   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())

plotbracket(0.1,0.13,0.32,0.025,'D')
plotbracket(0.23,0.255,0.28,0.025,'S')
plotbracket(0.18,0.41,0.25,0.025,'C')
plotbracket(0.67,0.79,0.25,0.025,'C')
plotbracket(0.67,0.79,0.25,0.025,'C')
plotbracket(0.27,0.73,0.215,0.025,'L')

plt.ylim(0.1,0.4)
# plt.xticks(list(np.round(np.arange(np.round(xmin,1),np.round(xmax,1)+1e-6,0.1),1)))

plt.yticks([])


plt.tight_layout()
plt.savefig('/home/tquah/Figures/slicealldata.png',dpi = 300)
