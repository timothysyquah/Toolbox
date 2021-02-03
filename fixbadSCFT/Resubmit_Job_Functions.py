#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 13:35:01 2021

@author: tquah
"""
import numpy as np
import os
import re
import pandas as pd
import itertools
from copy import deepcopy
from scipy.spatial.distance import cdist


def unique_category_values_F0dat(directory, file='F0_phases.dat'):
    data_path = os.path.join(directory, file)
    df = pd.read_csv(data_path, delimiter=' ', header=None)
    phases = list(df[0])
    plist = []
    for i in range(0, len(phases)):
        temp_tuple = (phases[i], float(df[1][i]), int(df[2][i]))
        plist.append(temp_tuple)

    category = directory.split('/')
    category_values = []
    for i in range(0, len(category)):
        extract_values = re.findall('[\d]*[.][\d]+', category[i])
        if len(extract_values) == 1:
            category_values.append(float(extract_values[0]))
        elif len(extract_values) > 1:
            t = ()
            for tup in extract_values:
                t = t + (float(tup),)
            category_values.append(t)
    category_values.append(plist)
    return category_values

def unique_category_values_diagnose(directory, file='diagnose.dat'):
    data_path = os.path.join(directory, file)
    df = pd.read_csv(data_path, delimiter=' ', header=None)
    phases = list(df[0])
    phaselist = []
    for i in range(0, len(phases)):
        temp_tuple = (phases[i], float(df[1][i]), int(df[2][i]), float(np.round(df[3][i], 3)))
        phaselist.append(temp_tuple)

    category = directory.split('/')
    category_values = []
    for i in range(0, len(category)):
        extract_values = re.findall('[\d]*[.][\d]+', category[i])
        if len(extract_values) == 1:
            category_values.append(float(extract_values[0]))
        elif len(extract_values) > 1:
            t = ()
            for tup in extract_values:
                t = t + (float(tup),)
            category_values.append(t)
    category_values.append(phaselist)
    return category_values

def Seperate_Categories_unique(listofcategory):
    depth = len(listofcategory[0]) - 2
    returnlist = []
    for i in range(0, depth):
        templist = []
        for j in range(0, len(listofcategory)):
            templist.append(listofcategory[j][i])

        templist = sorted(list(set(templist)))
        returnlist.append(templist)

    return returnlist

def prune_undesired_phases(lst,desired):
    newlist = []
    for node in lst:
        if node[-1] in desired:
            newlist.append(node)
    return newlist

def prune_index(distance,index,tol):
    loc = np.where(distance<tol)[0]
    return index[loc],distance[loc],loc

def Creating_Dictionary(dictionary, category, disflag=False,tolerance = 1e-6):
    key_tuple = ()
    header = ['Phase', 'FE', 'Status', 'Structure']
    for i in range(0, len(category) - 1):
        key_tuple += (category[i],)
    if disflag:
        header.append('DISFLAG')
        phases = []
        for i in range(0, len(category[-1]), 1):
            phases.append(category[-1][i][0])
        index_DIS = phases.index('DISPhase')
        for i in range(0, len(category[-1]), 1):
            free_energy_DIS = np.abs(category[-1][i][1] - category[-1][index_DIS][1])
            if free_energy_DIS < tolerance:
                category[-1][i] += (True,)
            else:
                category[-1][i] += (False,)
    dictionary[key_tuple] = pd.DataFrame.from_records(category[-1], columns=header)

def Combine_Dictionary(rawlist, disflag=False,tol = 1e-6):
    datadictionary = dict()
    for i in range(0, len(rawlist), 1):
        Creating_Dictionary(datadictionary, rawlist[i], disflag,tolerance=tol)
    return datadictionary

def Combine_List_Category(directorylist, unique_category_values=unique_category_values_diagnose):
    category_list = []
    for i in range(0, len(directorylist), 1):
        category_list.append(unique_category_values(directorylist[i]))
    return category_list

def DIS_Reporter(data_dictionary):
    DisList = []
    headers = list(data_dictionary)
    for i in range(0, len(headers)):
        DISFLAG_index = np.where(data_dictionary[headers[i]]['DISFLAG'] == True)
        list_of_DISPhases = list(data_dictionary[headers[i]]['Phase'][DISFLAG_index[0]])
        DIS_index = list_of_DISPhases.index('DISPhase')
        del list_of_DISPhases[DIS_index]
        if len(list_of_DISPhases) > 0:
            for j in range(0, len(list_of_DISPhases), 1):
                DisList.append(tuple((headers[i], list_of_DISPhases[j])))
    return DisList

def fA_Seperate_Bin(data_dictionary, ordered_category):
    res = list(itertools.product(*ordered_category[0:]))
    header = list(data_dictionary)
    fA_dict = dict()
    for i in range(0, len(header), 1):
        dictloc = res.index(list(data_dictionary)[i][0:-1])

        temp_header = list(fA_dict)

        if (res[dictloc] in temp_header) == False:
            fA_dict[res[dictloc]] = []
            fA_dict[res[dictloc]].append(header[i][-1])
        else:
            fA_dict[res[dictloc]].append(header[i][-1])

    for i in range(0, len(fA_dict), 1):
        fA_dict[list(fA_dict)[i]] = [np.array(sorted(list(set((fA_dict[list(fA_dict)[i]])))))]

    return fA_dict

def All_Phases(data_dictionary):
    unique_phases = []
    for i in data_dictionary:
        unique_phases += list(data_dictionary[i]['Phase'])
    return sorted(list(set(unique_phases)))

def prepfAdictionary(dictionary, phaselist):
    for i in dictionary:
        for phase in phaselist:
            fA_len = len(dictionary[i][0])
            dictionary[i].append(np.zeros(fA_len))


# Dictionary 2 is modified
def fA_populate(dictionary_1, dictionary_2, phaselist, flag):
    header1 = list(dictionary_1)
    header2 = list(dictionary_2)
    for i in range(0, len(header1)):
        index_fA_dict = header2.index(header1[i][0:-1])
        fA = header1[i][-1]
        fAlist = dictionary_2[header2[index_fA_dict]][0]
        # print(fAlist)
        index_innerarrays = np.where(fA == fAlist)[0][0]
        # print(index_innerarrays)
        assert len(phaselist) == len(dictionary_2[header2[index_fA_dict]]) - 1, 'Correct Number of Phases'
        phaselist_dictionary = list(dictionary_1[header1[i]]['Phase'])
        for j in range(1, len(phaselist) + 1):
            if phaselist[j - 1] in phaselist_dictionary:
                phase_index = phaselist_dictionary.index(phaselist[j - 1])
                dictionary_2[header2[index_fA_dict]][j][index_innerarrays] = dictionary_1[header1[i]][flag][phase_index]
            else:
                dictionary_2[header2[index_fA_dict]][j][index_innerarrays] = 4
    header = ['fA'] + phaselist
    for i in dictionary_2:
        dictionary_2[i] = pd.DataFrame(np.vstack(dictionary_2[i]).transpose(), columns=header)

def Setup_Grid(dictionary,rules):
    header = sorted(list(dictionary))
    corresponding_list = []
    array = 0
    for i in range(0, len(header), 1):
        temparray_data = np.array(dictionary[header[i]])
        shape = np.shape(temparray_data)
        base = np.ones((shape[0], len(header[0])))
        for j in range(0,len(header[0])):
            base[:, j] = rules[j](header[i][j])* base[:, j]
       
        temparray = np.hstack((base, temparray_data))
        if i == 0:
            array = deepcopy(temparray)
        else:
            array = np.vstack((array, temparray))
        for j in range(0, shape[0]):
            templist = list(header[i]) + [temparray_data[j, 0]]
            corresponding_list.append(templist)

    return array, corresponding_list

def ListtoArraytranslate(dislist, rules, phaselist):
    #    newlist = []
    disarray = np.zeros((len(dislist), len(dislist[0])))
    count = 0
    for i in dislist:
        index = phaselist.index(i[-1])
        disarray[count, 0] = rules[0](i[0]) 
        disarray[count, 1] = rules[1](i[1])
        disarray[count, 2] = i[2]
        disarray[count, 3] = index
        count += 1

    return disarray

# For DIS Logic
# 3 is ok, but it is an empty point
# 0 is good
# 1 is bad


# assume all criteria must be satisfied
def Nearest_Neighbor(x0_wphase, xgrid, criteria_arrays, criterias):
    x0 = x0_wphase[:, :-1]
    phases = x0_wphase[:, -1].astype(int)
    # assumes x0 is mx3 matrix and xgrid is a nx3 matrix
    nn_index = np.zeros_like(phases)
    nn_distance = np.zeros_like(phases).astype(float)
    for i in range(0, len(phases), 1):
        distance = cdist(xgrid,x0[i,:].reshape(1,len(x0[i,:])))
        new_order = np.argsort(distance.reshape(len(distance),))
        # print(new_order)
        search_logic = False
        count = 1
        while search_logic == False:
            index = new_order[count]
            logic_array = np.zeros(len(criteria_arrays))
            for j in range(0, len(criteria_arrays), 1):
                value = criteria_arrays[j][index, phases[i]]
                logic_array[j] = criterias[j](value)
            if np.product(logic_array) == False:
                count += 1
            else:
                search_logic = True
            if count == len(new_order):
                print('No Neighbors Found')
                break
        nn_index[i] = new_order[count]
        nn_distance[i] = distance[nn_index[i]]
    return nn_index, nn_distance

# this function is specific for my stuff, but can be changed by directory structure
def asym_name_format(somelist):
    return f'chiAB_{somelist[0]:0.4f}/NscA_{somelist[1][0]:0.1f}_NscB_{somelist[1][1]:0.1f}/fA{somelist[2]:0.5f}/{somelist[3]}'


def write_resubmit_file(file, list1, list2):
    file_open = open(file, 'w+')
    for i in range(0, len(list1)):
        text = list1[i] + ' ' + list2[i]
        file_open.write(text)
        file_open.write('\n')
    file_open.close()

def DIS_List_Gather(dictionary, phaselist):
    DIS_resubmitlist = []
    for i in dictionary:
        for phase in phaselist:
            testarray = np.array(dictionary[i][phase])
            # check point above
            logic_above = np.abs(testarray[2:] - testarray[1:-1])
            # check point below
            logic_below = np.abs(testarray[1:-1] - testarray[0:-2])
            combineLogic = logic_above + logic_below
            index = np.where((combineLogic == 1) & (testarray[1:-1] == 1))[0] + 1
            if len(index) > 0:
                for j in index:
                    templist = list(i) + [dictionary[i]['fA'][j]] + [phase]
                    DIS_resubmitlist.append(templist)
    return DIS_resubmitlist

def Divergent_List_Gather(dictionary, phaselist):
    header = list(dictionary)
    Div_list = []
    for i in range(0, len(dictionary), 1):
        subheader = list(dictionary[header[i]])[1:]
        fA = np.array(dictionary[header[i]]['fA'])
        array = np.array(dictionary[header[i]][phaselist])
        div_index = np.where(array == 1)
        if len(div_index[0]) > 0:
            for j in range(0, len(div_index[0])):
                templist = list(header[i]) + [fA[div_index[0][j]]] + [subheader[div_index[1][j]]]
                Div_list.append(templist)
    return Div_list

# def Structure_List_Gather(dictionary, plist, tol=0.9):
#     DIS_resubmitlist = []
#     for i in dictionary:
#         for phase in plist:
#             testarray = np.array(dictionary[i][phase])
#             logic_array = np.ones_like(testarray)
#             logicloc = np.where(testarray < tol)[0]
#             logic_array[logicloc] = 0
#             # check point above
#             logic_above = np.abs(logic_array[2:] - logic_array[1:-1])
#             # check point below
#             logic_below = np.abs(logic_array[1:-1] - logic_array[0:-2])
#             combineLogic = logic_above + logic_below
#             index = np.where((combineLogic == 1) & (testarray[1:-1] < tol))[0] + 1
#             if len(index) > 0:
#                 for j in index:
#                     templist = list(i) + [dictionary[i]['fA'][j]] + [phase]
#                     DIS_resubmitlist.append(templist)
#     return DIS_resubmitlist

def Structure_List_Gather(dictionary, plist, tol=0.9):
    DIS_resubmitlist = []
    for i in dictionary:
        for phase in plist:
            testarray = np.array(dictionary[i][phase])
            logic_array = np.ones_like(testarray)
            logicloc = np.where(testarray < tol)[0]
            index = deepcopy(logicloc)
            # logic_array[logicloc] = 0
            # # check point above
            # logic_above = np.abs(logic_array[2:] - logic_array[1:-1])
            # # check point below
            # logic_below = np.abs(logic_array[1:-1] - logic_array[0:-2])
            # combineLogic = logic_above + logic_below
            # index = np.where((combineLogic == 1) & (testarray[1:-1] < tol))[0] + 1
            if len(index) > 0:
                for j in index:
                    templist = list(i) + [dictionary[i]['fA'][j]] + [phase]
                    DIS_resubmitlist.append(templist)
    return DIS_resubmitlist

def phase_node_dictionary(dictionary,rules,depth):
    plotdict = dict()
    depth_loc = np.where(np.array(depth)>1)[0]
    if len(depth_loc)>1:
        print('Plotting not supported for dimensions higher than 2')
    else:
        for node in dictionary:
            nondis_phases = np.where(dictionary[node]['DISFLAG'] == False)[0]
            yvalue = rules[depth_loc[0]](node[depth_loc[0]])
            xvalue = node[-1]

            if len(nondis_phases) != 0:
                phase_loc = np.where(dictionary[node]['FE'][nondis_phases] ==
                                     np.min(dictionary[node]['FE'][nondis_phases]))[0]
                if len(phase_loc)>0:
                    
                    stable_phase = dictionary[node]['Phase'][nondis_phases[phase_loc[0]]]
                else:
                    stable_phase = 'DISPhase'

            else:
                stable_phase = 'DISPhase'

            if stable_phase not in list(plotdict):
                plotdict[stable_phase] = np.array([xvalue, yvalue])
            else:
                plotdict[stable_phase] = np.vstack((plotdict[stable_phase], np.array([xvalue, yvalue])))
    return plotdict, depth_loc

def special_node_dictionary(listofnodes,rules,depthloc):
    node_dict = dict()    
    for node in listofnodes:
        phase = node[-1]
        y = rules[depthloc[0]](node[depthloc[0]])
        x = node[-2]
        if phase not in list(node_dict):
            node_dict[phase] = np.array([x,y])
        else:
            node_dict[phase] = np.vstack((node_dict[phase],np.array([x,y])))
    return node_dict

def remove_duplicates(lst):
    tpls = [tuple(x) for x in lst]
    dct = list(dict.fromkeys(tpls))
    dup_free = [list(x) for x in dct]
    return dup_free

def prune_deep_dis(resub_array, dis_array, phaselist, rawlist, n=4):
    nphases = len(phaselist)
    dim = np.shape(resub_array)
    returnlist = [] 
    for i in range(0,dim[0]):
        distance = cdist(dis_array[:,0:nphases-1],\
                         resub_array[i,0:nphases-1].reshape((1,len(resub_array[i,0:nphases-1]))))
        order = np.argsort(distance.reshape(len(distance),))
        loc = order[1:n+1]
        if len(np.where(dis_array[loc,nphases:]==0)[0])>0:
            returnlist.append(rawlist[i])
    return returnlist
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
        



def Clean_Array(Array):
    unique_array = np.unique(Array[:,0:2],axis=0)
    # print(unique_array)
    shape = np.shape(unique_array)
    return_array = np.zeros((shape[0],3))
    for i in range(0,shape[0]):
        loc = np.where((unique_array[i,0]==Array[:,0])&(unique_array[i,1]==Array[:,1]))[0]
        if len(loc)==1:
            return_array[i,:] = Array[loc,:]
        if len(loc)>1:
            temp_compare = np.zeros(len(loc))
            count = 0
            for j in loc:
                temp_compare[count] = Array[j,2]
                count+=1
            choice = np.where(temp_compare==np.min(temp_compare))[0]
            if len(choice)>0:
                return_array[i,:] = Array[loc[choice[0]],:]
            else:
                return_array[i,:] = Array[loc[choice],:]
    return return_array

def Combine_Phases(dictionary,phase):  
    count = 0          
    header = list(dictionary)
    for key in header:
        if key.find(phase)!=-1:
#            print(key)
            if count == 0:
                temp_array = deepcopy(dictionary[key])
                count+=1
            else:
                temp_array = np.vstack((temp_array,dictionary[key]))
    return temp_array

