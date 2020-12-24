#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:35:01 2020

@author: tquah
"""



import numpy as np
import glob,re
import os,sys
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle
from scipy.interpolate import CubicSpline,interp1d,UnivariateSpline,barycentric_interpolate
from copy import deepcopy
import scipy as sp
import scipy.stats
import itertools

def data_extraction(path,filename,keyword):
    IDIR_og = os.getcwd()
    os.chdir(path)
    IDIR = os.getcwd()
    dir_dict = dict()
    for file in glob.glob(f'./chiAB_0.0289/ABratio_*/fA0.**/{filename}', recursive = True):
        dir_split = file.split('/')
        full_path = os.path.join(IDIR,file)
        op = open(full_path,'r')
        dataread = op.read().splitlines()
        op.close()
        dictionary_list = []
        for i in range(0,len(keyword)-1):
            dir1loc = dir_split.index([s for s in dir_split if keyword[i] in s][0])
            dir1=dir_split[dir1loc]
            dictionary_list.append(float(re.findall("\d+\.\d+", dir1)[0]))
    
        dir2loc = dir_split.index([s for s in dir_split if keyword[-1] in s][0])
        dir2=dir_split[dir2loc]
        innerstore = float(re.findall("\d+\.\d+", dir2)[0])
        
        for i in range(0,len(dataread),1):
            newlist = []
            splitdata = dataread[i].split(' ')
            phase = splitdata[0][:-5]
            newlist+=dictionary_list
            newlist.append(phase)      
            if int(splitdata[2])==2 or int(splitdata[2])==0 or int(splitdata[2])==3:
                if tuple(newlist) in dir_dict:
                    temparray = np.array([innerstore,float(splitdata[1])]) 
                    dir_dict[tuple(newlist)] = np.vstack((dir_dict[tuple(newlist)],temparray))
                else:
                    dir_dict[tuple(newlist)] = np.array([innerstore,float(splitdata[1])])
            else:
                print('Simulations have STATUS 0 1 3...revist extractF0.dat...')
                break
    os.chdir(IDIR_og)
    return dir_dict
def obtain_unique_dims(dictionary):
    
    keys = list(dictionary)
    var = len(keys[0])
    returnlist = []
    for i in range(0,var,1):
        temp = []

        for key in keys:
        
            temp.append(key[i])
        returnlist.append(sorted(list(set(temp))))
    return returnlist
        

path = "/media/tquah/TOSHIBA EXT/Projects/chiN_60_asymdir"
export_path = "/home/tquah/Figures/11-11-2020-PhasePlotterTool/"
filename = 'F0_phases.dat'
keyword = ['ABratio','f']
data_dict = data_extraction(path,filename,keyword)


header_sigma = list(data_dict)


for key in header_sigma:
    sort = np.argsort(data_dict[key][:,0])
    data_dict[key] = data_dict[key][sort]
    if key[1]=='LAM' or key[1]=='GYR' or key[1]=='FCC':
        data_dict.pop(key, None)
    if key[1].find('SIGMA')!=-1:
        if key[1]!='SIGMA':
            main_array = data_dict[(key[0],'SIGMA')]
            sigma_array = data_dict[key]
            compare_points = np.intersect1d(main_array[:,0],sigma_array[:,0],return_indices=True)
                
            if len(compare_points[0])>0:
                for i in range(0,len(compare_points),1):
                    if main_array[compare_points[1][i],1] > sigma_array[compare_points[2][i],1]:
                        main_array[compare_points[1][i],1] = sigma_array[compare_points[2][i],1]
            
            unique = np.setdiff1d(sigma_array[:,0],main_array[:,0])
            if len(unique)>0:
                indices = np.searchsorted(sigma_array[:,0], unique)
                main_array = np.vstack((main_array,sigma_array[indices,:]))
            sort = np.argsort(main_array[:,0])
            main_array = main_array[sort,:]
            data_dict[(key[0],'SIGMA')] = main_array
            data_dict.pop(key, None)

        
listofvar = obtain_unique_dims(data_dict)





rangemin = 0.10
rangemax = 0.50
c  = ['b','r','g','k','m','c']
markers=[',', '+', '*', '.', 'o', '*']
ratiomin = 1.0
ratiomax = 3.0
dF = 0.04
fArange = np.round(np.arange(rangemin,rangemax+1e-6,dF),2)
ABratio_range = np.arange(ratiomin,ratiomax+1e-6,0.5)
DIS_Bounds = [0.15,0.85]

#break from general methods
DISloc = listofvar[1].index('DIS')

plt.close('all')
data = []
phase_names = []

plt.figure()
Spline_Dictionary = dict()
for i in range(0,len(listofvar[0]),1):

    fA = data_dict[listofvar[0][i],'DIS'][:,0]

    H = data_dict[listofvar[0][i],'DIS'][:,1]
    indices = np.argsort(fA,0)
    Hsorted = H[indices]        
    fAsorted = fA[indices]
#    CS = CubicSpline(fAsorted,Hsorted)
    CS = interp1d(fAsorted,Hsorted,bounds_error=False)
    Hrange = CS(fArange)*1000
    if i==0:
        data_array = np.vstack((fArange,np.ones_like(fArange)*listofvar[0][i],Hrange)).transpose()
    if i!=0:
        temp_array = np.vstack((fArange,np.ones_like(fArange)*listofvar[0][i],Hrange)).transpose()
        data_array = np.vstack((data_array,temp_array))
        # plt.plot(temp_array)
    Spline_Dictionary[listofvar[0][i],'DIS'] = CS
    # plt.plot(fArange,Hrange)
    # plt.scatter(fAsorted,Hsorted)


# datalist = []
data = []
phase_names = []
for i in range(0,len(listofvar[1]),1):
    # plt.figure()
    if listofvar[1][i]=='LAM' or listofvar[1][i]=='GYR'or listofvar[1][i]=='FCC':
        continue
    # if i ==1:
    #     break
    # print(listofvar[1][i])
    
    if listofvar[1][i]=='DIS':
        data.append(data_array)
        phase_names.append(listofvar[1][i])

    else:
        count = 0
        if i!=0:
            del data_array_nondis
        for j in range(0,len(listofvar[0]),1):
            if (listofvar[0][j],listofvar[1][i]) in list(data_dict):
                
                # print((listofvar[0][j],listofvar[1][i]))
                dim = np.shape(data_dict[listofvar[0][j],listofvar[1][i]])
                
                if len(dim)>1:
                    fA = data_dict[listofvar[0][j],listofvar[1][i]][:,0]
                    H = data_dict[listofvar[0][j],listofvar[1][i]][:,1]

            
                    
                    indices = np.argsort(fA,0)
                    Hsort = H[indices]        
                    fAsort = fA[indices]
                    smallarealoc = np.where(fAsort<=0.6)[0]
                        
                        
                        

                    
                    
                    
                    if len(smallarealoc)>0:
                        fAsorted = fAsort[smallarealoc]
                        Hsorted = Hsort[smallarealoc]
                        
                        # print(x)
                        index = np.intersect1d(data_dict[listofvar[0][j],'DIS'][:,0],fAsorted,return_indices=True) 
                
                        yy = Hsorted-data_dict[listofvar[0][j],'DIS'][index[1],1]
                        loc = np.where(np.abs(yy)>1e-10)[0]
                        fA_actual = fAsorted[loc]
                        H_actual = Hsorted[loc]
#                        CS = CubicSpline(fAsorted,Hsorted)
                        CS = interp1d(fA_actual,H_actual,bounds_error=False)
                        Spline_Dictionary[listofvar[0][j],listofvar[1][i]] = CS
                        # CS = barycentric_interpolate(fA_actual,H_actual)

                        fA_local_minmax = np.array([np.min(fA_actual),np.max(fA_actual)])
                        local_loc = np.where((fArange >= fA_local_minmax[0]) & (fArange <= fA_local_minmax[1]))[0]
                        # print(fArange[local_loc])
                        
                        temp_array = np.zeros([len(fArange),3])
                        temp_array[:,0]=fArange
                        temp_array[:,1] = listofvar[0][j]
                        # print(local_fA)
                        temp_array[local_loc,2] = CS(fArange[local_loc])*1000
                        
                        
                        # if listofvar[0][j]=='BCC':
                        zeroloc = np.where(temp_array[:,2]==0)[0]
                        DISloc = np.where((temp_array[zeroloc,0]<DIS_Bounds[0]) | (temp_array[zeroloc,0]>DIS_Bounds[1]))

                        # temp_array[DISloc,2] = 1000*Spline_Dictionary[listofvar[0][j],'DIS'](temp_array[DISloc,0])+0.01

                        
                        
                        zeroloc = np.where(temp_array[:,2]==0)[0]

                        temp_array[zeroloc,2]=np.nan*fArange[zeroloc]
            else:
                temp_array = np.zeros([len(fArange),3])
                temp_array[:,0]=fArange
                temp_array[:,1] = listofvar[0][j]
                temp_array[:,2]=np.nan*fArange

                
                
            if count==0:
                data_array_nondis = temp_array
                count+=1
            elif count!=0:
                data_array_nondis = np.vstack((data_array_nondis,temp_array))


        data.append(data_array_nondis)
        phase_names.append(listofvar[1][i])
datascrub = deepcopy(data)

del data
data = []

for array in datascrub:
    # plt.figure()

    nan_array = np.isnan(array[:,2])
    not_nan_array = ~ nan_array
    # data.append(array[not_nan_array,:])
    dataphase = array[not_nan_array,:]
    fAunique = np.unique(dataphase[:,0])
    count2=0
    for i in range(0,len(fArange),1):
        fAloc = np.where(dataphase[:,0]==fArange[i])[0]
        
        if len(fAloc)>2:
            fA = dataphase[fAloc,0]
            ABratio = dataphase[fAloc,1]
            Hratio = dataphase[fAloc,2]
            
            CS = interp1d(ABratio,Hratio,bounds_error=False)
            ABratio_local_minmax = np.array([np.min(ABratio),np.max(ABratio)])
            local_loc = np.where((ABratio_range >= ABratio_local_minmax[0]) & (ABratio_range <= ABratio_local_minmax[1]))[0]
            ABratiolocal = ABratio_range[local_loc]
            Hlocal = CS(ABratiolocal)
            
            if count2 == 0:
                count2+=1
                data_array = np.ones((len(ABratio_range),3))*np.nan
                data_array[:,0] = fA[0]
                data_array[:,1] = ABratio_range
                data_array[local_loc,2] = Hlocal
                
                
                # data_array = np.vstack((np.ones_like(Hlocal)*fA[0],ABratiolocal,Hlocal)).transpose()
            else:
                # temp_array = np.vstack((np.ones_like(Hlocal)*fA[0],ABratiolocal,Hlocal)).transpose()
                temp_array = np.ones((len(ABratio_range),3))*np.nan
                temp_array[:,0] = fA[0]
                temp_array[:,1] = ABratio_range
                temp_array[local_loc,2] = Hlocal

                data_array = np.vstack((data_array,temp_array))
        else:
            if count2 == 0:
                count2+=1
                data_array = np.ones((len(ABratio_range),3))*np.nan
                data_array[:,0] = fA[0]
                data_array[:,1] = ABratio_range
                
                
                # data_array = np.vstack((np.ones_like(Hlocal)*fA[0],ABratiolocal,Hlocal)).transpose()
            else:
                # temp_array = np.vstack((np.ones_like(Hlocal)*fA[0],ABratiolocal,Hlocal)).transpose()
                temp_array = np.ones((len(ABratio_range),3))*np.nan
                temp_array[:,0] = fArange[i]
                temp_array[:,1] = ABratio_range
                data_array = np.vstack((data_array,temp_array))

    data.append(data_array)
    
    
    
    
    
    
    # plt.title(phase_names[count])
    # plt.scatter(data_array[:,0],data_array[:,1])
    # count +=1
# print(data)



# data = []
count = 0


    
for i in range(0,len(fArange)):
    for j in range(0,len(ABratio_range)):
        Emin = 50
        lowestphase = 'DIS'
        row = i*len(ABratio_range)+j
        # print(row)
        for k in range(0,len(data)):
            # print(data[k][i*len(ABratio_range)+j,2])
            if Emin>data[k][row,2]:
#                lowestphase = phase_names[k]
                Emin = data[k][row,2]
#                phaseindex = k
                
        for k in range(0,len(data)):
            if data[k][row,2]!=data[k][row,2]:
                data[k][row,2] = 0.05+(Emin)



# data = deepcopy(datascrub)
for i in range(0,len(data),1):
    # plt.figure()

    # x = data[i][:,0]
    # y = data[i][:,1]
    # plt.scatter(x,y)
    # plt.title(phase_names[i]+'Phase')
    # plt.xlim(0.09,0.51)
    # plt.ylim(0.9,3.1)
    
    
    name = phase_names[i]+'.dat'
    
    np.savetxt(name,data[i],fmt = '%0.5f',header = 'fA EPS H')



dim = int(np.ceil(np.sqrt(len(ABratio_range))))
fig, ax = plt.subplots(dim-1, dim)
k = 0
count = 0
llist = []
plist = []

for i in range(0,len(ABratio_range)):
    
    
    
    if i%dim==0:
        if count ==0:
            count+=1
        else:
            k+=1
    print(k,i%dim)
    
    ax[k,i%dim].set_title(('AB Ratio: %0.2f'%ABratio_range[i]))
    for j in range(0,len(data),1):
        if j!=2:
            
            xloc = np.where(data[j][:,1]==ABratio_range[i])[0]
            l, = ax[k,i%dim].plot(data[j][xloc,0],data[j][xloc,2]-data[2][xloc,2],label = phase_names[j])
            # ax[k,i%dim].set_ylim(-0.5,0.5)
            print(i)
            if i==0:
                llist.append(l)
                plist.append(phase_names[j])
    plt.legend(llist,plist)

plt.figure()


grid = np.zeros([len(fArange),len(ABratio_range)])
phase_count = np.zeros(len(phase_names))
for i in range(0,len(ABratio_range)):
    abratio = ABratio_range[i]
    for k in range(0,len(fArange),1):
        fA = fArange[k]
        Hstore = 1000
        phaseloc = 10
        for j in range(0,len(phase_names)):
#            print(j)
            data_array = data[j]
            fAlist = np.where(data_array[:,0]==fA)[0]
            ratiolist = np.where(data_array[:,1]==abratio)[0]
            index = np.intersect1d(fAlist, ratiolist)
            if Hstore>data_array[index,2]:
                Hstore = data_array[index,2]
                phase = phase_names
                phaseloc = j
        if phaseloc!=10:
            
            if phase_count[phaseloc]==0:
                plt.scatter(fA,abratio,color = c[phaseloc],marker = markers[phaseloc],s = 100.0,label = phase_names[phaseloc])
#                print(phase_names[phaseloc])
                phase_count[phaseloc]+=1
            else:
                plt.scatter(fA,abratio,color = c[phaseloc],marker = markers[phaseloc],s = 100.0)


plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.xlabel('fA')
plt.ylabel('$\epsilon$')

grid_figure = os.path.join(export_path,'grid.pdf')

plt.savefig(grid_figure,dpi = 300)