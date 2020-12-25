#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:46:23 2020

@author: tquah
"""


from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
import os
import re, glob
#paths to fields 1 and 2
#Fields_1 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields.in'
#Fields_2 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields_k.bin'


def my_function(X, Y): 
    return "Hello %s %s!" % (X, Y)

def Scaledown_Resolution(field,desired_npw):
    field_Re = field.AllFields.transpose()
    field_Im = field.AllFieldsImPart.transpose()
    coords = field.AllMeshCoords.transpose()
    npw = np.array(field.griddim)
    dim = len(npw)
    npw_ratio = npw/desired_npw
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
#    newfield_Im = field_Im[rowlist,:].transpose()
    return newfield_Re


def fields_compare(fields_1,fields_2,verbose = False):
    #Need to find lower reslution using shape and take every other field point to match grid!
    l1 = os.path.exists(fields_1)
    l2 = os.path.exists(fields_2)
    if l1==True and l2==True:
        fields1 = PolyFTSFieldReader()
        try:
            fields1.readFields(fields_1,verbose=verbose)
        except:
            print"*** Error during field reading ***"
        f1 = fields1.AllFields
        fields2 = PolyFTSFieldReader()
        try:
            fields2.readFields(fields_2,verbose=verbose)
        except:
            print "*** Error during field reading ***"
        f2 = fields2.AllFields
#        assert(np.shape(f1)[0]==np.shape(f2)[0])
#        assert(np.shape(f1)[1]==np.shape(f2)[1])
        if np.shape(f1)[1]>np.shape(f2)[1]:
            desired_npw = np.array(fields2.griddim)
            f1 = Scaledown_Resolution(fields1,desired_npw)
            
        elif np.shape(f2)[1]>np.shape(f1)[1]:
            desired_npw = np.array(fields1.griddim)
            f2 = Scaledown_Resolution(fields2,desired_npw)
        
        
        check_list = []
        for i in range(0,np.shape(f2)[0],1):
            f1norm1 = f1[i,:]/np.linalg.norm(f1[i,:])
            f2norm1 = f2[i,:]/np.linalg.norm(f2[i,:]) 
            check = np.dot(f1norm1,f2norm1)
            check_list.append(float(check))
        return check_list    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool to verify phases converged to the desired structure')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-l', '--LAMcount', action='store', default = 3,help='Number of LAM Checks')
    parser.add_argument('--ignorephase', action='store',default=[''], nargs='+',help='phases to ignore')
    parser.add_argument('--verbose', action='store',default=True,help='Verbose')

    parser.add_argument('-seedpath', '--seedpath', action='store', default='FIELD_CHECK',help='Field path 1')
    args = parser.parse_args()
#    args.fields1 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields.in'
#    args.fields2 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields_k.bin'
#    os.chdir("/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_asymBCC_fix/")
    WDIR = os.getcwd()
    fullseedpath = os.path.join(WDIR,args.seedpath)
    phases2ignore = [a+'Phase' for a in args.ignorephase ]
#    args.dirs = glob.glob("chiAB*/Nsc*/fA*")
    count = 0
    for directory in args.dirs:
#        print(directory)
        count+=1
        phaselist = glob.glob('*Phase')
        for phase in phaselist:
#            if phase=='HEXPhase':
#                so = open('STRUCTURE','w+')
#                so.write('2')
#                so.close()
#            else:
            if phase not in phases2ignore:
                phasefind = phase.find('Phase')
                phasename = phase[:phasefind]
                phase_directory = os.path.join(directory,phase)
                os.chdir(phase_directory)
                files = os.listdir('./')
                if 'fields_k.bin' in files:
    
                    if phasename=='LAM':
                        storevalue = np.zeros(args.LAMcount)
                        for i in range(0,args.LAMcount):                
                            
                            seedfield_name = phasename+str(i+1)+'_fields.in'
                            seedfield_path = os.path.join(fullseedpath,seedfield_name)
                            fieldoutput = fields_compare('fields_k.bin',seedfield_path)
                            storevalue[i] = np.min(np.abs(fieldoutput))
                        fieldvalue = np.max(storevalue)
                    else:
                        seedfield_name = phasename+'_fields.in'
                        seedfield_path = os.path.join(fullseedpath,seedfield_name)
                        fieldoutput = fields_compare('fields_k.bin',seedfield_path)
                        fieldvalue = np.min(np.abs(fieldoutput))
                    
                    so = open('STRUCTURE','w+')
                    if args.verbose:
                        print(directory)
                        print(phase+' %0.3f'%fieldvalue)
    
                    if fieldvalue>0.90:
                        so.write('2')
                    else:
                        so.write('0')
                    so.close()
                os.chdir(WDIR)
            
                        
                        
    
    
#    field_logic = fields_compare(args.fields1,args.fields2)

    
#print(fields_compare(Fields_1,Fields_2))