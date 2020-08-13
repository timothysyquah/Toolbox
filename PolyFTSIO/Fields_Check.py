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
#paths to fields 1 and 2
#Fields_1 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields.in'
#Fields_2 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields_k.bin'


def my_function(X, Y): 
    return "Hello %s %s!" % (X, Y)



def fields_compare(fields_1,fields_2):
    l1 = os.path.exists(fields_1)
    l2 = os.path.exists(fields_2)
    if l1==True and l2==True:
        fields1 = PolyFTSFieldReader()
        try:
            fields1.readFields(fields_1,True)
        except:
            print"*** Error during field reading ***"
        f1 = fields1.AllFields
        fields2 = PolyFTSFieldReader()
        try:
            fields2.readFields(fields_2,True)
        except:
            print "*** Error during field reading ***"
        f2 = fields2.AllFields
        assert(np.shape(f1)[0]==np.shape(f2)[0])
        assert(np.shape(f1)[1]==np.shape(f2)[1])
        check_list = []
        for i in range(0,np.shape(f2)[0],1):
            f1norm1 = f1[i,:]/np.linalg.norm(f1[i,:])
            f2norm1 = f2[i,:]/np.linalg.norm(f2[i,:]) 
            check = np.dot(f1norm1,f2norm1)
            check_list.append(float(check))
        return check_list    


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool to verify phases converged to the desired structure')
    parser.add_argument('-f1', '--fields1', action='store', default='',help='Field path 1')
    parser.add_argument('-f2', '--fields2', action='store', default='',help='Field path 2')
    args = parser.parse_args()
#    args.fields1 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields.in'
#    args.fields2 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields_k.bin'
    field_logic = fields_compare(args.fields1,args.fields2)

    
#print(fields_compare(Fields_1,Fields_2))