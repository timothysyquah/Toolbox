#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 00:37:47 2021

@author: tquah
inputfile=__INFILE__
outputfile=__OUTFILE__
polyftsdir=__PFTPATH__

"""


import numpy as np
import re 
import glob 
import os
import math
import shutil
from copy import deepcopy
import subprocess


def round_up_to_even(f):
    return math.ceil(f / 2.) * 2

def OpenGetData(file):
    op = open(file,'r')
    data = op.read()
    op.close()
    return data


def Checkstatus(fullist,phase):
    IDIR = os.getcwd()
    dellist = []
    for i in range(0,len(fullist)):
        fullpath = os.path.join(fullist[i],phase)
        os.chdir(fullpath)
        status = int(OpenGetData('STATUS'))
        if status!=2:
            index = fullist.index(fullist[i])
            dellist.append(index)
        os.chdir(IDIR)
    # print(dellist)
    for i in range(len(dellist)):
        fullist.pop(dellist[-i])

def ListPattern(string,list_):
    r = re.compile(string)
    newlist = list(filter(r.match,list_))
    return newlist

def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 

           
def PurgeList(fulllist,Nsc,Nbblist):    
    list1 = []
    nbblist = []
    for i in range(0,len(fulllist),1):
        splitlist = fulllist[i].split('/')
        ncslist = ListPattern("nsc",splitlist)
        nbblist = ListPattern("nbb",splitlist)
        Nscvalue = float(extract_value(ncslist[0])[0])
        Nbbvalue = float(extract_value(nbblist[0])[0])
        if Nsc==Nscvalue: #if nsc is desired
            if Nbbvalue+1 in Nbblist: #if nbb is desired
                list1.append(fulllist[i])
    return list1
    #compile nbb
def delete_elements(inlist,string):
    delete_list = []
    for i in range(len(inlist)):
        if inlist[i]==string:
            delete_list.append(i)
    delete_list.reverse()
    for i in range(len(delete_list)):
        inlist.pop(delete_list[i])
    

def ParseInfile(path,infile,listofproperties):
    IDIR = os.getcwd()
    os.chdir(path)
    data = OpenGetData(infile)
    datalist = data.splitlines()
    datadict = dict()
    for i in range(len(listofproperties)):
        datadict[listofproperties[i]] = []
    
    
    for i in range(len(datalist)):
        for j in range(len(listofproperties)):
            if listofproperties[j] in datalist[i]:
                datadict[listofproperties[j]].append(extract_value(datalist[i]))
    os.chdir(IDIR)
    return datadict


def ExtractDomain(path,outfile):
    IDIR = os.getcwd()
    os.chdir(path)
    data = OpenGetData(outfile)
    datasplit = data.splitlines()
    for i in range(len(datasplit)):
        if "Final simulation cell:" in datasplit[i]:
            values = extract_value(datasplit[i])
            D_0 = float(values[1])*10**float(values[2])
    os.chdir(IDIR)
    return D_0


def Custom_Parser(d_dictionary,ReplaceList):
    NewDicitionary = dict()
    NewDicitionary['GaussSmearWidth'] = [float(d_dictionary['GaussSmearWidth'][0][0]),ReplaceList[0]]
    NewDicitionary['Nbb'] = [float(d_dictionary['nbeads'][0][0]),ReplaceList[1]]
    NewDicitionary['NbbA'] = [float(d_dictionary['nperblock'][0][0]),ReplaceList[2]]
    NewDicitionary['NbbB'] = [float(d_dictionary['nperblock'][0][1]),ReplaceList[3]]
    NewDicitionary['Nsc'] = [float(d_dictionary['nbeads'][1][0]),ReplaceList[4]]
    NewDicitionary['Narms'] = [float(d_dictionary['numarms'][0][0]),ReplaceList[5]]
    NewDicitionary['Astart'] = [float(d_dictionary['backbonegraftingstart'][0][0]),ReplaceList[6]]
    NewDicitionary['Aend'] = [float(d_dictionary['backbonegraftingend'][0][0]),ReplaceList[7]]
    NewDicitionary['Bstart'] = [float(d_dictionary['backbonegraftingstart'][1][0]),ReplaceList[8]]
    NewDicitionary['Bend'] = [float(d_dictionary['backbonegraftingend'][1][0]),ReplaceList[9]]
    NewDicitionary['Lx'] = [float(d_dictionary['Lx'][0]),ReplaceList[10]]
    NewDicitionary['npwx'] = [float(d_dictionary['Lx'][0])/NewDicitionary['GaussSmearWidth'][0],ReplaceList[12]]
    NewDicitionary['npwyz'] = [float(d_dictionary['npwyz'][0]),ReplaceList[13]]
    NewDicitionary['Lyz'] = [float(d_dictionary['npwyz'][0])*NewDicitionary['GaussSmearWidth'][0],ReplaceList[11]]
    NewDicitionary['chiN'] = [float(d_dictionary['chiN12'][0][1]),ReplaceList[14]]
    NewDicitionary['zetaN'] = [float(d_dictionary['compressibility_invzetaN'][0][0]),ReplaceList[15]]
    NewDicitionary['C'] = [float(d_dictionary['C'][0]),ReplaceList[16]]
    return NewDicitionary
    
def make_submit(SUBMITFILE,PHASE,JOBNAME,WDIR,PFTSPATH):
    with open('%s/submit.sh' % WDIR, 'w') as fout:
        with open(SUBMITFILE,'r') as f:
            for line in f:
                line = line.replace('__JOBNAME__',JOBNAME)
                line = line.replace('__INFILE__',f'{PHASE}.in')
                line = line.replace('__PFTPATH__',PFTSPATH)
                line = line.replace('__OUTFILE__',f'{PHASE}.out')
                fout.write(line)
    fout.close
    f.close()
    
def editfile(outfile,infile,dictionary,WDIR):
    header = list(dictionary)    
    with open(f'%s/{outfile}' % WDIR, 'w') as fout:
        with open(infile,'r') as f:
            for line in f:
                for i in range(len(header)):
                    line = line.replace(dictionary[header[i]][1],str(dictionary[header[i]][0]))
                fout.write(line)
    fout.close
    f.close()
def Convert_List_To_Path(list_,absolute=True):
    if absolute:
        path = ''
        for i in range(0,len(list_),1):
            path +=list_[i]+'/'
    else:
        path = ''
        for i in range(1,len(list_),1):
            path +=list_[i]+'/'
    return path

def Create_New_Path(path,string,string_replace):
    segment_path = path.split('/')
    newlist = [string_replace if x==string else x for x in segment_path]
    newpath = Convert_List_To_Path(newlist,absolute=True)
    return newpath



def MakeDirectoryStructure(wdir):
    logicrun = 1    
    try:
    # Create target Directory
        os.makedirs(wdir)
        # print("Directory " , wdir ,  " Created ") 
        logicrun = 1
    except FileExistsError:
        # print("Directory " , wdir ,  " already exists")        
        logicrun = 0
    return logicrun



ReplaceList = ['__a__',\
               '__nbb__','__nbbA__','__nbbB__',\
               '__nsc__','__narms__',\
                   '__Astart__','___Aend__',\
               '__Bstart__','__Bend__',\
                '__Lx__','__Lyz__','__npwx__','__npwyz__',\
                    '__chiN__','__zeta__','__c__']
        
ParameterList = ['GaussSmearWidth',\
               'nbeads','nperblock','numarms',\
                   'backbonegraftingstart','backbonegraftingend',\
                    'chiN12','compressibility_invzetaN']
        

    
IDIR = os.getcwd()

#################################################################################    
#Parameters
#################################################################################    
#parse old lam scripts
npw_yz = 32
C = 2.0
Nsc = 20.0
Nbb = [50, 70, 100, 130, 170, 210]
low = 4
high = 10
dL = 2
submit_jobs = False
# Nsc = 40.0
# Nbb = [50, 60, 80, 100, 120, 140]


#################################################################################    
#Pathing
#################################################################################
Directory_Structure = 'SCFT/f*/chi*/nsc*/nbb*'
field_name = 'fields_k.bin'
SubmitPath = 'SEEDS/submitGPU.sh'
phase = 'LAM3D'
LAMPath = f'SEEDS/{phase}_CL_template.in'
polyfts_path = '/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release'
#################################################################################
absolute_main_path = os.path.join(IDIR,Directory_Structure)
fulldirectory_list = glob.glob(absolute_main_path)
Checkstatus(fulldirectory_list,f'{phase}Phase')
infile = f"{phase}.in"
list_of_jobs = PurgeList(fulldirectory_list,Nsc,Nbb)

for i in range(len(list_of_jobs)):
    path = os.path.join(list_of_jobs[i],f"{phase}Phase") 
    field_path = os.path.join(path,field_name)
    data_dictionary = ParseInfile(path,f'{phase}.in',ParameterList)
    D0 = ExtractDomain(path,f'{phase}.out')
    D0_round = round_up_to_even(D0)
    data_dictionary['Lx'] = [D0_round]
    data_dictionary['C'] = [C]
    data_dictionary['npwyz'] = [npw_yz]
    final_dictionary = Custom_Parser(data_dictionary,ReplaceList)
    Dlist = np.arange(D0_round-low,D0_round+high+1e-6,dL,dtype =int)
    # npwxlist = Dlist/final_dictionary['GaussSmearWidth'][0]
    newpath = Create_New_Path(path,'SCFT','CL')
    
    for j in range(len(Dlist)):
        print(newpath)

        wdir = os.path.join(newpath,f'L_{Dlist[j]}/')
        logic = MakeDirectoryStructure(wdir)
        
        if logic==0:
            print('Skipping!-Already Exists')
            continue  
        submit_full_path = os.path.join(IDIR,SubmitPath)
        make_submit(submit_full_path,phase,'CL-BB',wdir,polyfts_path)
        lam_full_path = os.path.join(IDIR,LAMPath)
        final_dictionary['Lx'][0] = Dlist[j]
        final_dictionary['npwx'][0] = Dlist[j]/final_dictionary['GaussSmearWidth'][0]
        editfile(f'{phase}.in',lam_full_path,final_dictionary,wdir)
        field_outpath = os.path.join(wdir,'fields.in')
        shutil.copyfile(field_path,field_outpath)
        if submit_jobs:
            os.chdir(wdir)
            cmd="qsub submit.sh"
            #cmd="sbatch submit.sh"
            subprocess.call(cmd.split())
            os.chdir(IDIR)
        
        
        
