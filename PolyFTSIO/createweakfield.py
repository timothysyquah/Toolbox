#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 00:49:50 2021

@author: tquah
"""
import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/')
import iotools as io
import numpy as np
from copy import deepcopy
from PolyFTSIO import *
def getfield_stats(field):
    return np.mean(field),np.std(field)

# ##Method 1 use random noise to perturb structure
# # Read input file
# infile = '/home/tquah/Projects/CL_SEEDS/LAM_field.bin'
# outfile = '/home/tquah/Projects/CL_SEEDS/testlam.dat'
# coords, fields = io.ReadBinFile(infile)
# field_list = []
# for i in range(len(fields[0,0,0,:])):
#     if i%2==0:
        
#         length =len(np.ravel(fields[:,:,:,i]))
#         mean,std =getfield_stats(fields[:,:,:,i])
#         random_vars = np.random.normal(0,std,length)
        
#         reshape_random = \
#             np.reshape(random_vars,\
#                        (np.int(np.round(length**(1/3))),\
#                         np.int(np.round(length**(1/3))),\
#                             np.int(np.round(length**(1/3)))))
#         field_new = fields[:,:,:,i]+reshape_random
#         field_list.append(field_new)
#     # else:
#     #     field_list.append(field_new)

# fieldoutput = np.stack(field_list,axis = 3)    

# io.WriteDatFile(outfile, coords, fieldoutput, iskspace = False, iscomplex = False)

# Method 2 take average of two fields one DIS one LAM

infilelist = ['/home/tquah/Projects/CL_SEEDS/averagedis.bin','/home/tquah/Projects/CL_SEEDS/LAM_field.bin']
outfile = '/home/tquah/Projects/CL_SEEDS/avglamwithavgdis.dat'

weight = [30 ,1] # x:1

# coords, fields = io.ReadBinFile(infile)
fields = []
coords = []
for i in range(len(infilelist)):
    coord_local, fields_local = io.ReadBinFile(infilelist[i])
    fields.append(fields_local*weight[i])
    coords.append(coord_local)
    
newfield = np.stack(fields)
avgfield = np.sum(fields,axis = 0)/np.sum(weight)
#get real fields 
fieldout = avgfield[:,:,:,[0,2]]


io.WriteDatFile(outfile, coord_local, fieldout, iskspace = False, iscomplex = False)

# Method 3 lets take an average of 100 random fields then average it with 

# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/Average_Random_Field'
# n = 100
# outfile = '/home/tquah/Projects/CL_SEEDS/avgdis.dat'

# for i in range(1,n+1):
#     print(i)
#     newpath = os.path.join(path,f'n_{i}/fields.bin')
#     coord_local, fields_local = io.ReadBinFile(newpath)

#     if i==1:
#         coordmain = deepcopy(coord_local)
#         fieldmain = deepcopy(fields_local)
#     else:
#         fieldmain+=fields_local
        
        
# io.WriteDatFile(outfile, coordmain, fieldmain/n, iskspace = False, iscomplex = False)






