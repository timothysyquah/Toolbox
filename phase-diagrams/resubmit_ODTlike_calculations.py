#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 23:06:09 2020

@author: tquah
"""

import numpy as np
import os 
import shutil
import re
import glob
import time
# import joblib
import pandas as pd
import itertools
from copy import deepcopy

def unique_category_values(directory,file='diagnose.dat'):
    data_path = os.path.join(directory,file)
    df = pd.read_csv(data_path,delimiter=' ',header=None)
    phases = list(df[0])
    phaselist = []
    for i in range(0,len(phases)):
        temp_tuple = (phases[i],float(df[1][i]),int(df[2][i]))
        phaselist.append(temp_tuple)
        
    category = directory.split('/')
    category_values = []
    for i in range(0,len(category)):
        extract_values = re.findall('[\d]*[.][\d]+' ,category[i])
        if len(extract_values)==1:
            category_values.append(float(extract_values[0]))
        elif len(extract_values)>1:
            t = ()
            for tup in extract_values:
                t = t + (float(tup),)
            category_values.append(t)
    category_values.append(phaselist)
    return category_values

    
    
# print(category)
    # for i in range(0,depth,1):
    

def Seperate_Categories_unique(listofcategory):
    depth = len(listofcategory[0])-2
    returnlist = []
    for i in range(0,depth):
        templist = []
        for j in range(0,len(listofcategory)):
           templist.append(listofcategory[j][i])
        
        templist = sorted(list(set(templist)))
        returnlist.append(templist)
        
    return returnlist
    

def Creating_Dictionary(dictionary,category,disflag = False):
    key_tuple = ()
    header = ['Phase', 'FE', 'Status']
    for i in range(0,len(category)-1):
        key_tuple+=(category[i],)
    if disflag:
        header.append('DISFLAG')
        phases = []
        for i in range(0,len(category[-1]),1):
            phases.append(category[-1][i][0])
        index_DIS = phases.index('DISPhase')
        for i in range(0,len(category[-1]),1):
            free_energy_DIS = np.abs(category[-1][i][1]-category[-1][index_DIS][1])
            if free_energy_DIS<1e-8:
                category[-1][i]+=(True,)
            else:
                category[-1][i]+=(False,)
    dictionary[key_tuple] = pd.DataFrame.from_records(category[-1], columns =header) 

    # return dictionary

def Combine_Dictionary(rawlist,disflag = False):
    datadictionary = dict()
    for i in range(0,len(rawlist),1):
        Creating_Dictionary(datadictionary,rawlist[i],disflag)
    return datadictionary
def Combine_List_Category(directorylist):
        category_list = []
        for i in range(0,len(directorylist),1):
            category_list.append(unique_category_values(directorylist[i]))
        return category_list

def DIS_Reporter(data_dictionary):
    DisList = []
    headers = list(data_dictionary)
    for i in range(0,len(headers)):
        DISFLAG_index = np.where(data_dictionary[headers[i]]['DISFLAG']==True)
        list_of_DISPhases = list(data_dictionary[headers[i]]['Phase'][DISFLAG_index[0]])
        DIS_index = list_of_DISPhases.index('DISPhase')
        del list_of_DISPhases[DIS_index]
        if len(list_of_DISPhases)>0:
            for j in range(0,len(list_of_DISPhases),1):
                DisList.append(tuple((headers[i],list_of_DISPhases[j])))
    return DisList

def fA_Seperate_Bin(data_dictionary,ordered_category):
    res = list(itertools.product(*ordered_category[0:])) 
    header = list(data_dictionary)
    fA_dict = dict()
    for i in range(0,len(header),1):
        dictloc = res.index(list(data_dictionary)[i][0:-1])
        
        temp_header = list(fA_dict)
        
        if (res[dictloc] in temp_header) == False:
            fA_dict[res[dictloc]] = []
            fA_dict[res[dictloc]].append(header[i][-1])
        else:
            fA_dict[res[dictloc]].append(header[i][-1])
    
    for i in range(0,len(fA_dict),1):
        fA_dict[list(fA_dict)[i]] = [np.array(sorted(list(set(( fA_dict[list(fA_dict)[i]])))))]
        
    return fA_dict
def All_Phases(data_dictionary):
    unique_phases = []
    for i in data_dictionary:
        unique_phases+=list(data_dictionary[i]['Phase'])
    return sorted(list(set(unique_phases)))

def prepfAdictionary(dictionary,phaselist):
    for i in dictionary:
        for phase in phaselist:
            fA_len = len(dictionary[i][0])
            dictionary[i].append(np.zeros(fA_len))
        

#Dictionary 2 is modified
def fA_populate(dictionary_1,dictionary_2,phaselist,flag):
    header1 = list(dictionary_1)
    header2 = list(dictionary_2)
    for i in range(0,len(header1)):
        index_fA_dict = header2.index(header1[i][0:-1])
        fA = header1[i][-1]
        fAlist = dictionary_2[header2[index_fA_dict]][0]
        # print(fAlist)
        index_innerarrays = np.where(fA==fAlist)[0][0]
        # print(index_innerarrays)
        assert len(phaselist)==len(dictionary_2[header2[index_fA_dict]])-1,'Correct Number of Phases'
        phaselist_dictionary =list(dictionary_1[header1[i]]['Phase'])
        for j in range(1,len(phaselist)+1):
            if phaselist[j-1] in phaselist_dictionary:
                phase_index = phaselist_dictionary.index(phaselist[j-1])
                dictionary_2[header2[index_fA_dict]][j][index_innerarrays] = dictionary_1[header1[i]][flag][phase_index]
            else:
                dictionary_2[header2[index_fA_dict]][j][index_innerarrays] = 4
    header = ['fA']+phaselist
    for i in dictionary_2:
        dictionary_2[i] = pd.DataFrame(np.vstack(dictionary_2[i]).transpose(), columns =header) 


def DIS_List_Gather(dictionary,phaselist):
    DIS_resubmitlist = []
    for i in dictionary:
        for phase in phaselist:
            testarray = np.array(dictionary[i][phase])
            #check point above
            logic_above = np.abs(testarray[2:]-testarray[1:-1])
            #check point below
            logic_below =  np.abs(testarray[1:-1]-testarray[0:-2])
            
            combineLogic = logic_above+logic_below
            index = np.where((combineLogic==1)&(testarray[1:-1]==1))[0]+1
            if len(index)>0:
                for j in index:
                    templist = list(i)+[dictionary[i]['fA'][j]]+[phase]
                    DIS_resubmitlist.append(templist)
                    
    return DIS_resubmitlist

def Setup_Grid(dictionary)     :  

    header = sorted(list(dictionary))
    corresponding_list = []
    for i in range(0,len(header),1):
        temparray_data = np.array(dictionary[header[i]])
        shape = np.shape(temparray_data)
        base = np.ones((shape[0],len(header[0])))
        
        base[:,0] = header[i][0]*2100*base[:,0]
        base[:,1] = np.sqrt((header[i][1][1]+1)/(1+header[i][1][0]))*base[:,1]
        
        temparray = np.hstack((base,temparray_data))
        
        if i==0:
            array = deepcopy(temparray)
        else:
            array = np.vstack((array,temparray))    
        for j in range(0,shape[0]):
                
            
            templist =  list(header[i])+[temparray_data[j,0]]
            corresponding_list.append(templist)
            
    return array,corresponding_list

def ListtoArraytranslate(dislist,phaselist):
#    newlist = []
    disarray = np.zeros((len(dislist),len(dislist[0])))    
    count = 0
    for i in dislist:
        
        index = phaselist.index(i[-1])
        disarray[count,0] = i[0]*2100
        disarray[count,1] = np.sqrt((i[1][1]+1)/(i[1][0]+1))
        disarray[count,2] = i[2]
        disarray[count,3] = index
        count+=1
    
    return disarray
#For DIS Logic
#3 is ok, but it is an empty point
#0 is good
#1 is bad


def DISLogic(value):
    if value ==4 or value ==1:
        return False
    elif value ==0:
        return True
    else:
        print('Warning Wrong Criteria')
        return False
    
#For Status Logic
# 1 is bad
# 0/3 is ok (depending)
# 2 is good


def StatusLogic_converged(value):
    if value ==3 or value ==0 or value==1 or value==4:
        return False
    elif value ==2:
        return True
    else:
        print('Warning Wrong Criteria')
        return False

def StatusLogic(value):
    if value ==3 or value ==0 or value==2:
        return True
    elif value ==1 or value==4:
        return False
    else:
        print('Warning Wrong Criteria')
        return False

#assume all criteria must be satisfied
def Nearest_Neighbor(x0_wphase,xgrid,criteria_arrays,criterias):
    x0 = x0_wphase[:,:-1]
    phases = x0_wphase[:,-1].astype(int)
    #assumes x0 is mx3 matrix and xgrid is a nx3 matrix
    nn_index = np.zeros_like(phases)
    nn_distance = np.zeros_like(phases).astype(float)
    for i in range(0,len(phases),1):
        xgrid_distance = np.abs(xgrid-x0[i,:])   
        #use eucliden distance 
        distance = np.sqrt(np.dot(x0[i,:],xgrid_distance.transpose()))
        new_order = np.argsort(distance)
        # print(new_order)
        search_logic = False
        count = 1
        while search_logic==False:
            index = new_order[count]
            logic_array = np.zeros(len(criteria_arrays))
            for j in range(0,len(criteria_arrays),1):
                value = criteria_arrays[j][index,phases[i]]
                logic_array[j] = criterias[j](value)
            if np.product(logic_array)==False:
                count+=1
            else:
                search_logic=True
            if count ==len(new_order):
                print('No Neighbors Found')
                break
        nn_index[i] = new_order[count]
        nn_distance[i] = distance[nn_index[i]]
    return nn_index,nn_distance

#this function is specific for my stuff, but can be changed by directory structure
def asym_name_format(somelist):
    return f'chiAB_{somelist[0]:0.4f}/NscA_{somelist[1][0]:0.1f}_NscB_{somelist[1][1]:0.1f}/fA{somelist[2]:0.5f}/{somelist[3]}'
   
def write_resubmit_file(file,list1,list2):
    file_open = open(file,'w+')
    for i in range(0,len(list1)):
        text = list1[i]+' '+list2[i]
        file_open.write(text)
        file_open.write('\n')
    file_open.close()
    
def Divergent_List_Gather(dictionary,phaselist):
    header = list(dictionary)
    Div_list = []
    for i in range(0,len(dictionary),1):
        subheader = list(dictionary[header[i]])[1:]
        fA = np.array(dictionary[header[i]]['fA'])
        array = np.array(dictionary[header[i]][phaselist])
        div_index = np.where(array==1)
        if len(div_index[0])>0:
            for j in range(0,len(div_index[0])):
                templist = list(header[i])+[fA[div_index[0][j]]]+[subheader[div_index[1][j]]]
                Div_list.append(templist)
    return Div_list



#Path to Some Directory
os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_asymBCC_fix/')

wdir = os.getcwd()
dirs = glob.glob("chiAB_*/Nsc*/fA*")



#do not include DIS
phaselist = ['BCCPhase','LAMPhase','HEXPhase','GYRPhase']

# Want a general way to place data in a searchable dictionary that can be used for other scripts




#determine folders depth importance

#first determine depth from first and last directory

depth = len(dirs[0].split('/'))
assert len(dirs[-1].split('/'))==depth,'Inconsistant Directory Structure'

#load data into python
listofcategory = Combine_List_Category(dirs)
ordered_category = Seperate_Categories_unique(listofcategory)
#organize dictionary
data_dictionary = Combine_Dictionary(listofcategory,disflag = True)
#create fA lists for each category 
unique_phases = All_Phases(data_dictionary)


##This is for DIS Flag
#determine if you want the phases specified or all phases found 
fA_dict = fA_Seperate_Bin(data_dictionary,ordered_category)
# Grid it out
prepfAdictionary(fA_dict,phaselist)
fA_populate(data_dictionary,fA_dict,phaselist,'DISFLAG')
##This is for status flag
fA_status = fA_Seperate_Bin(data_dictionary,ordered_category)
# Grid it out
prepfAdictionary(fA_status,phaselist)
fA_populate(data_dictionary,fA_status,phaselist,'Status')
###Note If I want to do additional flags I need to modify combine dictionary function and create an additional flag like DIS
#Example in the future I would like a structure check which would show up here!
#Example in the future I would also like a STATUS Existance Check which would show up here!


#Here it checks for DIS
DIS_resubmitlist = DIS_List_Gather(fA_dict,phaselist)
Divergent_resubmitlist = Divergent_List_Gather(fA_status,phaselist)
combined_resubmitlist = DIS_resubmitlist+Divergent_resubmitlist

#Here it checks for Status
allarray,corresponding_list = Setup_Grid(fA_dict)
statusarray,corresponding_list2 = Setup_Grid(fA_status)


#Here it creates Grids
DISsarray = ListtoArraytranslate(combined_resubmitlist,phaselist)  
gridarray = allarray[:,0:len(phaselist)-1]
criteria_array = [allarray[:,len(phaselist)-1:],statusarray[:,len(phaselist)-1:]]


#Note here we can attach any number of crteria to check and grid out
criterion = [DISLogic,StatusLogic]
distance = Nearest_Neighbor(DISsarray,gridarray,criteria_array,criterion)

#Here we find the nearest neighbor to the point
nn_index,nn_distance = Nearest_Neighbor(DISsarray,gridarray,criteria_array,criterion)
seed_list = []

#We are just making the listt of seeds
count = 0
for i in nn_index:
    seed_list.append(corresponding_list[i]+[combined_resubmitlist[count][-1]]) 
    count+=1

#prepare real names and writing a file that shows the seed directory and the directory in question:
# We use a bash script from here to resubmit jobs and move fields/change box guesses
resubmitList_names = []
SeedList_names = []
for i in range(0,len(seed_list)):
    resubmitList_names.append(asym_name_format(combined_resubmitlist[i]))
    SeedList_names.append(asym_name_format(seed_list[i]))
write_resubmit_file('Resubmit_DIS.dat',resubmitList_names,SeedList_names)
dataread = open('Resubmit_DIS.dat','r')
testvalidate = dataread.read()
dataread.close()
    
    