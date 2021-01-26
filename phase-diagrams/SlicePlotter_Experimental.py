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
plt.figure(figsize=(7,5))
yspace = np.linspace(0,0.29,100)
data_array = np.loadtxt('asymslice.dat')
for i in range(0,len(data_array[:,0])):
    plt.plot(data_array[i,0]*np.ones_like(yspace),yspace,'--k',alpha = 0.3)
plt.yticks(color='w')
plt.xlabel('$f_A$')
overallphases = ['LAM','CYL','SPH']

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
symbol = ['x','+','+','x']

color = ['k','y','g','b']
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
        print(phase)
        refindex= [i for i, x in enumerate(fullreference) if x ==ref]
        phaseindex= [i for i, x in enumerate(fullphases) if x ==phase]
        height = overallphases.index(phase)
        loc = list(set(refindex) & set(phaseindex))
        if len(loc)>0:
            
            # print(Nsc_Backbone[loc])
            
            

            mloc = np.where(Nsc_Backbone[loc]<1.0)[0]
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
                    plt.scatter(xfilter,(height/10+0.05)*np.ones_like(xfilter),c=color[count1],marker = symbol[count2],label = name)      
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
arraylist  = [np.array([0.241]),\
              np.array([0.3,0.781])]
label = ['SPH','CYL']
symbol = ['o','x','x','^']

i=0

for array in arraylist:
    height = overallphases.index(label[i])

    plt.scatter(array,(height/10+0.05)*np.ones_like(array),c= 'r',marker = symbol[i],label = label[i]+'-Present Work')
    i+=1
plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

xmin = 0.08
xmax = 0.7
xarray = np.ones(6)
xarray[0] = xmin
sortdata = np.sort(data_array[:,0])
xarray[1:-1] = sortdata[:4]
xarray[-1] = xmax
textsize = 20
order = ['D','S','C','G','L']
avg = lambda x: (x[1]+x[0])/2

y= 0.3
for i in range(0,len(xarray)-1,1):
    xpos = avg(xarray[i:i+2])
    plt.text(xpos,y,order[i],color = 'k',size = textsize,ha = 'center')



# plt.xlim(xmin,xmax)




plt.ylim(0.0,0.35)
# plt.xticks(list(np.round(np.arange(np.round(xmin,1),np.round(xmax,1)+1e-6,0.1),1)))

plt.yticks([])
plt.tight_layout()
plt.savefig('/home/tquah/Figures/slice_alicedata.png',dpi = 300)
