#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  9 10:16:35 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import os
pathlist = ['/home/tquah/Projects/positions/LAM/0.0104/density_chain0.dat',\
            '/home/tquah/Projects/positions/LAM/0.0189/density_chain0.dat',\
            '/home/tquah/Projects/positions/LAM/0.0289/density_chain0.dat',\
            '/home/tquah/Projects/positions/LAM/0.0389/density_chain0.dat']
                        
export_path = '/home/tquah/Projects/positions/LAM/Analysis'
plt.close('all')
errortotal = dict()


# fig, axs = plt.subplots(2, , sharex=True, sharey=True)

pathcount = 0

for path in pathlist:
    export_name = path.split('/')[-2][-4:]+'.pdf'
    export_name_raw = path.split('/')[-2][-4:]+'raw.pdf'
    export_name_14 = path.split('/')[-2][-4:]+'14.pdf'

    data = np.loadtxt(path)
    # fig = plt.figure()
    # ax = plt.subplot(111)
    
    shape = np.shape(data)
    # count =0
    zero_1 = np.max(data[:,1])
    zero_2 = np.max(data[:,2])
    errortotal[float(path.split('/')[-2])]= []
    # for i in range(0,shape[1],1):
    #     # print(i)
    #     if np.max(data[:,i])<zero_1:
    #         count+=1
    #         # print(count)
    
    #         if count%4==0:
    #             # print(i)
    #             diff = (data[:,i]-data[:,1])**2
    #             errortotal[float(path.split('/')[-2])].append(np.sum(diff))
    #             ax.plot(data[:,0],diff,label = ('%d'%(count)))
    
    # ax.legend(bbox_to_anchor=(1.1, 1.05))
    # # plt.colorbar()
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    # plt.xlabel('Position (x)')
    # plt.ylabel(r'$ \left( \rho_{bb}-\rho_{sc,i} \right )^2 $')
    # plt.tight_layout()
    
    # pathsave = os.path.join(export_path,export_name)
    # plt.savefig(pathsave,dpi=300)
    fig = plt.figure()
    ax = plt.subplot(111)
    
    

    count =0
    halfmax = 0
    peak = []
    for i in range(0,shape[1],1):
        # print(i)
        if np.max(data[:,i])<zero_1:
            count+=1
            # print(count)
    
            if count<=12 and count%2:
                # print(i)  
                diff = (data[:,i]-data[:,1])**2
                ax.plot(data[:,0],data[:,i],label = ('%d'%(count)))
                peak.append(np.max(data[:,i]))
                middle_loc = np.where(np.max(data[:,i])==data[:,i])[0][0]
                middle = data[middle_loc,0]
                # print(middle)
                if np.max(data[:,i])>halfmax:
                    halfmax = np.max(data[:,i])
    
    ax.legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.ylim(np.min(peak)*0.99,np.max(peak)*1.01)
    plt.xlabel('Position (x)')
    plt.ylabel(r'$\rho$')
    plt.title(('$\chi = %0.4f$'%float(path.split('/')[-2])))
    plt.xlim(middle*0.85,middle*1.15)
    plt.tight_layout()

    pathsave = os.path.join(export_path,export_name_14)
    plt.savefig(pathsave,dpi=300)

    fig = plt.figure()
    ax = plt.subplot(111)
    
    shape = np.shape(data)
    count =0
    zero_1 = np.max(data[:,1])
    zero_2 = np.max(data[:,2])
    
    for i in range(0,shape[1],1):
        # print(i)
        if i==0:
            ax.plot(data[:,0],data[:,1],label = 'Backbone')
        # if i
        if pathcount==0 and i==20:
            ax.plot(data[:,0],data[:,i],label = ('%d'%(20)))

        if np.max(data[:,i])<=zero_1:
            count+=1
            # print(count)
            # print(i)
            

            if count%4==0:
                print(i)
                ax.plot(data[:,0],data[:,i],label = ('%d'%(count)))
                if np.max(data[:,i])>halfmax:
                    halfmax = np.max(data[:,i])
    
    plt.ylim(halfmax/2,halfmax*1.1)
    plt.xlim(np.max(data[:,0])/2,np.max(data[:,0]))

    ax.legend()
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.xlabel('Position (x)')
    plt.ylabel(r'$\rho$')
    plt.title(('$\chi = %0.4f$'%float(path.split('/')[-2])))
    
    plt.tight_layout()
    
    pathsave = os.path.join(export_path,export_name_raw)
    plt.savefig(pathsave,dpi=300)
    pathcount+=1
    # break