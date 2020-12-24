#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 14:12:35 2020

@author: tquah
"""


# op = open('asymslice.dat','r')
# data = op.read()
# op.close()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.close('all')
plt.figure()
yspace = np.linspace(0,1,100)
data_array = np.loadtxt('asymslice.dat')
for i in range(0,len(data_array[:,0])):
    plt.plot(data_array[i,0]*np.ones_like(yspace),yspace,'k')
plt.yticks(color='w')
plt.xlabel('$f_A$')


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
symbol = ['h','+','+','h']
color = ['k','r','g','b']
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
        refindex= [i for i, x in enumerate(fullreference) if x ==ref]
        phaseindex= [i for i, x in enumerate(fullphases) if x ==phase]
        
        
        loc = list(set(refindex) & set(phaseindex))
        if len(loc)>0:
            
            # print(Nsc_Backbone[loc])
            
            

            mloc = np.where(Nsc_Backbone[loc]<0.8)[0]
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
                
                print(ref)
                print(floc)

                
                
                # print(yfilter)
                # nscnbbloc = np.where(Nsc_Backbone[loc]<1.0)[0]
                # print(nscnbbloc)
            
                # print(xfilter)
            
                name = phase+'-'+ref
                if len(xfilter)>0:
                    plt.scatter(xfilter,0.5*np.ones_like(xfilter),c=color[count1],marker = symbol[count2],label = name)      
        count2+=1
    count1+=1

ypos = 6
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.xlim(0,0.5)
plt.savefig('/home/tquah/Figures/slice_otherdata.pdf',dpi = 300)

plt.figure()

for i in range(0,len(data_array[:,0])):
    plt.plot(data_array[i,0]*np.ones_like(yspace),yspace,'k')
plt.yticks(color='w')
plt.xlabel('$f_A$')
arraylist  = [np.array([0.384,0.457]),\
               np.array([0.241]),\
                   np.array([0.3,1-0.651,1-0.716,1-0.771]),\
                      np.array([1-0.883])]
label = ['LAM','BCC','HEX','DIS']
symbol = ['+','o','h','^']

i=0

for array in arraylist:
    plt.scatter(array,0.5*np.ones_like(array),marker = symbol[i],label = label[i])
    i+=1
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
plt.xlim(0,0.5)

plt.savefig('/home/tquah/Figures/slice_alicedata.pdf',dpi = 300)
