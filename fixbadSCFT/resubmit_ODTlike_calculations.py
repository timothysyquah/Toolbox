#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 23:06:09 2020

@author: tquah
"""
import os 
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
toolpath = op.read()
op.close()
chainpath =toolpath+'/phase-diagrams/'
import numpy as np
import sys
sys.path.append(chainpath)
import os
import glob
# matplotlib.use('TKAgg')
import matplotlib.pyplot as plt
import datetime
from Resubmit_Job_Functions import * 
import matplotlib.colors as colors
from copy import deepcopy
##Import Colors and Marker List


def DISLogic(value):
    if value == 4 or value == 1:
        return False
    elif value == 0:
        return True
    else:
        print('Warning Wrong Criteria')
        return False

# For Status Logic
# 1 is bad
# 0/3 is ok (depending)
# 2 is good


def StatusLogic_converged(value):
    if value == 3 or value == 0 or value == 1 or value == 4:
        return False
    elif value == 2:
        return True
    else:
        print('Warning Wrong Criteria')
        return False

def StatusLogic(value):
    if value == 3 or value == 0 or value == 2:
        return True
    elif value == 1 or value == 4:
        return False
    else:
        print('Warning Wrong Criteria')
        return False

def StructureLogic(value,tol = 0.9):
    if value >= tol:
        return True
    else:
        return False


colors_list = list(colors._colors_full_map.values())
marker = ['+', 'o', 's','v','^','<','>','x','p','h','H','*']
#plt.close('all')
## Path to Some Directory
#os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_asymBCC_fix/')
os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_asymBCC_fix/')
exportpath = '/home/tquah/Figures/DOCUMENTATION/'
phasediagram_tag = 'chiN_fA'
phasetolerance = 0.90
# os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected')
# os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected')
os.chdir('/media/tquah/TOSHIBA EXT/Projects/sweep-asym-armlength_BCC_fix/')

#Keywords you need to use
wdir = os.getcwd()
dirs = glob.glob("chiAB*/Nsc*/fA*")
output_name = 'RESUBMIT.dat'
# do not include DIS
desired_phaselist = ['BCCPhase', 'HEXPhase','LAMPhase','GYRPhase']
# desired_phaselist = ['GYRPhase']
desired_phaselist = ['HEXPhase']
# desired_phaselist = ['LAMPhase']
desired_phaselist = ['BCCPhase', 'HEXPhase','LAMPhase']

# desired_phaselist = ['BCCPhase']

visualize_submissions = True

#Rules-some need to be tweaked if you do not want seeding via y dimension penalize a rule harder!
tuplerule = lambda x: np.sqrt((x[1]+1)/(x[0]+1))
chiNrule = lambda x : x*20
rules = [chiNrule,tuplerule]

# determine folders depth importance
# first determine depth from first and last directory
depth = len(dirs[0].split('/'))
assert len(dirs[-1].split('/')) == depth, 'Inconsistant Directory Structure'

# load data into python
listofcategory = Combine_List_Category(dirs)
ordered_category = Seperate_Categories_unique(listofcategory)
depth_category = [len(x) for x in ordered_category]
# organize dictionary
data_dictionary = Combine_Dictionary(listofcategory, disflag=True,tol = 1e-6)
# create fA lists for each categoryarray
unique_phases = All_Phases(data_dictionary)
phaselist = deepcopy(unique_phases)
phaselist.remove('DISPhase')
##This is for DIS Flag
# determine if you want the phases specified or all phases found
fA_dict = fA_Seperate_Bin(data_dictionary, ordered_category)
# Grid it out
prepfAdictionary(fA_dict, phaselist)
fA_populate(data_dictionary, fA_dict, phaselist, 'DISFLAG')
##This is for status flag
fA_status = fA_Seperate_Bin(data_dictionary, ordered_category)
# Grid it out
prepfAdictionary(fA_status, phaselist)
fA_populate(data_dictionary, fA_status, phaselist, 'Status')
## This is for structure flag
fA_structure = fA_Seperate_Bin(data_dictionary, ordered_category)
# Grid it out
prepfAdictionary(fA_structure, phaselist)
fA_populate(data_dictionary, fA_structure, phaselist, 'Structure')


# Collects Points that fail the following criteria
DIS_resubmitlist = DIS_List_Gather(fA_dict, phaselist) #only detects if at least one surrounding point is NOT DIS
Divergent_resubmitlist = Divergent_List_Gather(fA_status, phaselist) #detects 1, but can be modified for 0 and 3
Structure_resubmitlist = Structure_List_Gather(fA_structure, phaselist,tol =phasetolerance) #detects with a tolerance of 0.9, but uses similar logic to DIS

# combine all the things I want to resubmit
# combined_resubmitlist = DIS_resubmitlist + Divergent_resubmitlist + Structure_resubmitlist
combined_resubmitlist = DIS_resubmitlist #+ Divergent_resubmitlist + Structure_resubmitlist

#remove duplicates
# combined_resubmitlist = remove_duplicates(combined_resubmitlist)
# combined_resubmitlist = prune_undesired_phases(combined_resubmitlist,desired_phaselist)
# Here it checks for Status
allarray, corresponding_list = Setup_Grid(fA_dict,rules)
statusarray, corresponding_list2 = Setup_Grid(fA_status,rules)
structurearray, corresponding_list3 = Setup_Grid(fA_structure,rules)


# Here it creates Grids
Resubmitarray = ListtoArraytranslate(combined_resubmitlist,rules, phaselist)

#Prunes points surrounded by DIS
combined_resubmitlist = prune_deep_dis(Resubmitarray,allarray,\
                                       phaselist,combined_resubmitlist)
Resubmitarray = ListtoArraytranslate(combined_resubmitlist,rules, phaselist)

gridarray = allarray[:, 0:len(phaselist) - 1]
full_criteria_array = [allarray,statusarray,structurearray]
criteria_array = [allarray[:, len(phaselist) - 1:], statusarray[:, len(phaselist) - 1:],
                  structurearray[:, len(phaselist) - 1:]]

# Note here we can attach any number of crteria to check and grid out
criterion = [DISLogic, StatusLogic, StructureLogic]
# distance = Nearest_Neighbor(Resubmitarray, gridarray, criteria_array, criterion)

# Here we find the nearest neighbor to the point
nn_index, nn_distance = Nearest_Neighbor(Resubmitarray, gridarray, criteria_array, criterion)


def prune_index(distance,index,tol):
    loc = np.where(distance<tol)[0]
    return index[loc],distance[loc],loc

nn_index,nn_distance,nn_loc = prune_index(nn_distance,nn_index,0.03)


final_prune = []
for i in nn_loc:
    final_prune.append(combined_resubmitlist[i])
del combined_resubmitlist
combined_resubmitlist = deepcopy(final_prune)



seed_list = []

# We are just making the listt of seeds
count = 0
for i in nn_index:
    seed_list.append(corresponding_list[i] + [combined_resubmitlist[count][-1]])
    count += 1

# prepare real names and writing a file that shows the seed directory and the directory in question:
# We use a bash script from here to resubmit jobs and move fields/change box guesses
resubmitList_names = []
SeedList_names = []
for i in range(0, len(seed_list)):
    resubmitList_names.append(asym_name_format(combined_resubmitlist[i]))
    SeedList_names.append(asym_name_format(seed_list[i]))
write_resubmit_file(output_name, resubmitList_names, SeedList_names)
dataread = open(output_name, 'r')
testvalidate = dataread.read()
dataread.close()

#Visualize what is going on!
if visualize_submissions:
    plotdict, depth_loc = phase_node_dictionary(data_dictionary,rules,depth_category)
    resubmit_dict = special_node_dictionary(combined_resubmitlist,rules,depth_loc)
    seed_dict = special_node_dictionary(seed_list,rules,depth_loc)
    phasekey = list(set(list(plotdict)+list(resubmit_dict)+list(seed_dict)))
    for phase in phasekey:
        fig, ax = plt.subplots(1, figsize = (5,5), sharex='all')
        count = 0
        for phase_background in plotdict:
            ax.scatter(plotdict[phase_background][:, 0], \
                       plotdict[phase_background][:, 1],\
                       marker=marker[count], color = colors_list[count],\
                       label = phase_background, alpha = 0.5, s = 60.0,zorder=1)
            count += 1
        if phase in resubmit_dict:    
            ax.scatter(resubmit_dict[phase][:, 0], resubmit_dict[phase][:, 1],\
               marker='o', color = 'r',s = 10.0,zorder=3,label = 'Resubmit')
            ax.scatter(seed_dict[phase][:, 0], seed_dict[phase][:, 1],\
                       marker='o', color = 'k', s = 10.0,zorder=3,label='Seed')
            for i in range(0,len(resubmit_dict[phase]),1):
                xpair = np.array([resubmit_dict[phase][i, 0],seed_dict[phase][i, 0]])
                ypair = np.array([resubmit_dict[phase][i, 1],seed_dict[phase][i, 1]])
                ax.plot(xpair,ypair,'r')
        ax.set_title(phase)
        ax.legend(bbox_to_anchor=(1.01, 0.5), loc="center left", borderaxespad=0)
        plt.tight_layout()
        today = datetime.datetime.now()
        export_name = f'{today:%Y-%m-%d}_{phase}_{phasediagram_tag}.png'
        full_path = os.path.join(exportpath,export_name)
        plt.savefig(full_path,dpi = 300)
