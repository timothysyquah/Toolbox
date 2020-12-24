#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 08:30:17 2020

@author: tquah
"""
from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
import os
import re, glob
import matplotlib.pyplot as plt
from copy import deepcopy
#from scipy.interpolate import 
#writePolyFTSDatFile(outfilename, griddim, hcell, complexdata, kspacedata, fielddata, fielddataimpart=[])
import_path = '/home/tquah/TestDIR/readwrite_different_resolutions/HEX_fields.in'
export_path = '/home/tquah/TestDIR/readwrite_different_resolutions/fields.out'
field.readFields(import_path,verbose=True)
field_Re = field.AllFields.transpose()
field_Im = field.AllFieldsImPart.transpose()
coords = field.AllMeshCoords.transpose()
npw = np.array(field.griddim)

#
desirednpw = np.array([16,16])
dim = len(npw)
npw_ratio = npw/desirednpw
#
coordsgrid = [] #np.zeros_like(coords)

newcoords = []
for i in range(0,np.shape(coords)[1]):
    coordshape = np.reshape(coords[:,i],tuple(npw))
    coordsgrid.append(coordshape)
newcoordsgrid = []
if dim==1:
    newcoordsgrid.append(coordsgrid[0][::npw_ratio[0]])
elif dim==2:
    newcoordsgrid.append(coordsgrid[0][::npw_ratio[0],::npw_ratio[1]])
    newcoordsgrid.append(coordsgrid[1][::npw_ratio[0],::npw_ratio[1]])
elif dim==3:
    newcoordsgrid.append(coordsgrid[0][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
    newcoordsgrid.append(coordsgrid[1][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
    newcoordsgrid.append(coordsgrid[2][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
else:
    print "Dimensions Not supported"

count = 0
for gridlist in newcoordsgrid:
    
    
    if count==0:
        smallcoords = np.ravel(gridlist)
    else:
        smallcoords = np.vstack((smallcoords,np.ravel(gridlist)))
    count+=1
smallcoords = smallcoords.transpose()
rowlist = []
for i in range(0,np.shape(smallcoords)[0],1):
    if dim==1:
        row = smallcoords[i]
    else:
        row = tuple(smallcoords[i,:])

    loc = np.where((coords==row).all(axis=1))[0][0]
    rowlist.append(loc)


newfield_Re = field_Re[rowlist,:].transpose()
newfield_Im = field_Im[rowlist,:].transpose()

def Scaledown_Resolution(field,desired_npw):
    field_Re = field.AllFields.transpose()
    field_Im = field.AllFieldsImPart.transpose()
    coords = field.AllMeshCoords.transpose()
    npw = np.array(field.griddim)
    dim = len(npw)
    npw_ratio = npw/desirednpw
    coordsgrid = []
    for i in range(0,np.shape(coords)[1]):
        coordshape = np.reshape(coords[:,i],tuple(npw))
        coordsgrid.append(coordshape)
    newcoordsgrid = []
    if dim==1:
        newcoordsgrid.append(coordsgrid[0][::npw_ratio[0]])
    elif dim==2:
        newcoordsgrid.append(coordsgrid[0][::npw_ratio[0],::npw_ratio[1]])
        newcoordsgrid.append(coordsgrid[1][::npw_ratio[0],::npw_ratio[1]])
    elif dim==3:
        newcoordsgrid.append(coordsgrid[0][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
        newcoordsgrid.append(coordsgrid[1][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
        newcoordsgrid.append(coordsgrid[2][::npw_ratio[0],::npw_ratio[1],::npw_ratio[2]])
    else:
        print "Dimensions Not supported"
    count = 0
    for gridlist in newcoordsgrid:
        if count==0:
            smallcoords = np.ravel(gridlist)
        else:
            smallcoords = np.vstack((smallcoords,np.ravel(gridlist)))
        count+=1
    smallcoords = smallcoords.transpose()
    rowlist = []
    for i in range(0,np.shape(smallcoords)[0],1):
        if dim==1:
            row = smallcoords[i]
        else:
            row = tuple(smallcoords[i,:])
    
        loc = np.where((coords==row).all(axis=1))[0][0]
        rowlist.append(loc)
    newfield_Re = field_Re[rowlist,:].transpose()
    newfield_Im = field_Im[rowlist,:].transpose()
    return newfield_Re,newfield_Im
