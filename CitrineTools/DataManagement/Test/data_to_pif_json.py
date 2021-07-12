#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:34:42 2021

@author: tquah
"""

import numpy as np
import os
import re
import math as m
import time
import networkx as nx
import matplotlib.pyplot as plt
import string
plt.close('all')
from pysmiles import write_smiles
from pypif import pif
from pypif.obj import *


def nCr(n,r):
    return int(m.factorial(n)/m.factorial(r)/(m.factorial(n-r)))

class chain():
    def __init__(self,statistic,n,species_order,composition):
        self.n = n
        self.species_order = species_order
        self.composition = composition
        self.statistic = statistic
    #returns properties in a list
    def return_properties(self):
        all_members = self.__dict__.keys()
        
        
        return [ (item, self.__dict__[item]) for item in all_members]

class sidechain(chain):
    def __init__(self,statistic,n,species_order,composition,graft_start,graft_end,graft_number):
        super().__init__(statistic,n,species_order,composition)
        self.graft_start =  graft_start
        self.graft_end = graft_end
        self.graft_number = graft_number
        

def OutputParser(text):
    DataDictionary = dict()
    for line in text:
        if 'initialized in' in line:
            d = int(re.findall(r'\d+', line)[0])
            break
    arch_list = []
    spacegroup_name = 'None'

    for i in range(len(text)):
        if 'Space group name' in text[i]: 
            spacegroup_name = text[i][text[i].index('=')+1:].strip()
        if 'Real-space collocation mesh' in text[i]: 
            npwarray = np.zeros(d)
            npw = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i])
            for k in range(d):
                npwarray[k] = int(npw[k])

        if "Number of monomer species" in text[i]:
            specs = int(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i])[0])

        if 'Chi parameters (original):' in text[i]: 
            combo = nCr(specs,2)
            chimatrix = np.zeros([combo,combo])
            for j in range(combo):
                values = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+j+1])
                for k in range(combo):
                    chimatrix[j,k] = float(values[k])
                    
        
        if 'Final simulation cell' in text[i]:
            celltensor = np.zeros([d,d])
            for j in range(d):
                values = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+j])
                for k in range(d):
                    a = float(values[2*k+1])
                    b = float(values[2*k+2])
                    finval = a*10**b
                    celltensor[j,k] = finval
            hamval = float(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i-2])[0])
        
        
        
        if 'Discrete Backbone' in text[i]:
            stat = text[i+1]
            loc = stat.index('+')
            stat = stat[loc+1:].strip()
            
            speclist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+4])
            specarray = np.zeros(len(speclist),dtype = int)
            for j in range(len(speclist)):
                specarray[j] = int(speclist[j])
            
            blist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+5])
            barray = np.zeros(len(blist))
            for j in range(len(blist)):
                barray[j] = float(blist[j])
            beadlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+6])
            beadarray = np.zeros(len(beadlist),dtype = int)
            for j in range(len(beadlist)):
                beadarray[j] = int(beadlist[j])

            mainback = [stat,specarray,barray,beadarray]
            arch_list.append(mainback)
            
            
            
        if 'Side arm index' in text[i]:
            stat = text[i+1]
            loc = stat.index('+')
            stat = stat[loc+1:].strip()
            
            speclist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+4])
            specarray = np.zeros(len(speclist),dtype = int)
            for j in range(len(speclist)):
                specarray[j] = int(speclist[j])
            
            blist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+5])
            barray = np.zeros(len(blist))
            for j in range(len(blist)):
                barray[j] = float(blist[j])
            beadlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+6])
            beadarray = np.zeros(len(beadlist),dtype = int)
            for j in range(len(beadlist)):
                beadarray[j] = int(beadlist[j])
            graftlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+9])
            graftarray = np.zeros(len(graftlist),dtype = int)
            for j in range(len(graftlist)):
                graftarray[j] = int(graftlist[j])

            mainback = [stat,specarray,barray,beadarray,graftarray]
            arch_list.append(mainback)

    DataDictionary['Space_Group'] = spacegroup_name
    DataDictionary['NPW'] = npwarray
    DataDictionary['CHI'] = chimatrix
    DataDictionary['ARCH'] = arch_list
    DataDictionary['CELL'] = celltensor
    DataDictionary['FREE_ENERGY'] = hamval
    return DataDictionary




#lets start with 
# A-A-A-A-A-B-B-B-B-B-B
# | | | | | | | | | | |
# A A A A A B B B B B B
# | | | | | | | | | | |
# A A A A A B B B B B B
# | | | | | 
# A A A A A 
# | | | | | 
# A A A A A 


        
def node_edge(n,comp,species,basenode=0,graft=False,graftpoint=None):
        n_array = np.arange(1,n+1e-6,1,dtype=int)+basenode
        ends = n*np.array(comp)
        
        start = 0
        end = 0
        edges = np.zeros([n-1,2])
        edges[:,0] = n_array[0:-1]
        edges[:,1] = n_array[1:]
        edgelist = list(map(tuple, edges))
        nodelist = []
        for i in range(0,len(comp),1):
            end += ends[i]
            nodelist.append(list(n_array[int(start):int(end)]))
            start+=ends[i]
        if graft:
            edgelist.append(tuple((graftpoint,n_array[0])))
        return nodelist,edgelist
class molecule():
    def __init__(self,listofchains):
        self.sidechain = []
        for chain in listofchains:
            
            if type(chain).__name__=='chain':
                self.backbone = chain
            elif type(chain).__name__=='sidechain': 
                self.sidechain.append(chain)
            else:
                print('Warning Chain Object Not Supported')
    def create_network(self,plot=True):
        #first create backbone
        label = list(string.ascii_lowercase)

        nbb = vars(self.backbone)['n']
        composition = vars(self.backbone)['composition']
        species= vars(self.backbone)['species_order']
        nodelist,edgelist = node_edge(nbb,composition,species)
        G = nx.Graph()
        for i in range(0,len(composition),1):
            G.add_nodes_from(nodelist[i],element = label[species[i]-1])
        G.add_edges_from(edgelist)
        # second create sidechains 
        node_count = nbb
        if len(self.sidechain)>0:
            for i in range(0,len(self.sidechain),1):
                nsc = vars(self.sidechain[i])['n']
                composition_sc = vars(self.sidechain[i])['composition']
                species_sc= vars(self.sidechain[i])['species_order']
                graftnumber = vars(self.sidechain[i])['graft_number']
                graftstart = vars(self.sidechain[i])['graft_start']
                graftend = vars(self.sidechain[i])['graft_end']                
                graft_array = np.linspace(graftstart,graftend,graftnumber,dtype = int)
                for i in range(0,len(graft_array),1):                
                    nodesclist,edgesclist = node_edge(nsc,composition_sc,\
                                                      species_sc,node_count,graft=True,\
                                                          graftpoint=graft_array[i])      
        
                    for j in range(0,len(composition_sc),1):
                        G.add_nodes_from(nodesclist[j],element = label[species_sc[j]-1])
                        node_count += len(nodesclist[j])

                    G.add_edges_from(edgesclist)

        
        # if plot:    
        #     nx.draw(G, with_labels=True, font_weight='bold')
        return G



# chainlist.append(chain('DGC',20,[1,2,3],[0.2,0.6,0.2]))
# chainlist.append(sidechain('DGC',3,[1],[1.0],1,4,4))
# chainlist.append(sidechain('DGC',3,[2],[1.0],5,16,12))
# chainlist.append(sidechain('DGC',3,[3],[1.0],17,20,4))

def list_convert_object(lst):
    chainlist = []
    for i in range(len(lst)):
        if i==0:
            Ntot = np.sum(lst[i][3]) 
            chainlist.append(chain(lst[i][0],Ntot,lst[i][1],lst[i][3]/Ntot))
        else:
            chainlist.append(sidechain(lst[i][0],lst[i][3][0],\
                                       [lst[i][1][0]],[1.0],int(np.min(lst[i][4]))+1,\
                                           int(np.max(lst[i][4]))+1,len(lst[i][4])))
    return chainlist

def dict_to_pif(Dictionary):

    chemical_system = ChemicalSystem()
    #START to Create Objects
    chemical_system.chemical_formula = 'N_{sc}='+str(Dictionary['NSC'])+'-f_{A/C}='+str(Dictionary['fAfC'])+'_'+Dictionary['PHASE']
    
    SMILES = Property()
    SMILES.name = 'SMILES'
    SMILES.scaler = Dictionary['SMILES']
    SMILES.datatype="Computational"

    
    RESOLUTION = Property()
    RESOLUTION.name = 'Number of Planewaves'
    RESOLUTION.vectors =  Dictionary['NPW'].tolist()
    RESOLUTION.datatype="Computational"
    
    NSC = Property()
    NSC.name = 'Side Chain Length'
    if len(Dictionary['ARCH'])==1:
       NSC.scaler = 0
    else:
        NSC.scaler = Dictionary['ARCH'][1][3][0]
    NSC.datatype="Computational"
  
    
    CS = Property()
    CS.name = 'Chain Statistics'
    CS.scaler = Dictionary['ARCH'][0][0]
    CS.datatype="Computational"

    SPACEGroup = Property()
    SPACEGroup.name = 'Space Group'
    SPACEGroup.scaler = Dictionary['Space_Group']
    SPACEGroup.datatype="Computational"

    FreeEnergy = Property()
    FreeEnergy.name = 'Free Energy'
    FreeEnergy.scaler = Dictionary['FREE_ENERGY']
    FreeEnergy.datatype="Computational"
 
    B = Property()
    B.name = 'Statistical Segment Length'
    B.vectors = Dictionary['ARCH'][0][2].tolist()
    B.datatype="Computational"
    B.units = 'b/sqrt(6)'
    CHI = Property()
    CHI.name = 'Flory-Huggin Interaction Matrix'
    CHI.matrices = Dictionary['CHI'].tolist()
    CHI.datatype="Computational"

    CELL = Property()
    CELL.name = 'Cell Tensor'
    CELL.matrices = Dictionary['CELL'].tolist()
    CELL.units = 'b/sqrt(6)'
    CELL.datatype="Computational"
    fAfC = Property()
    fAfC.name = 'Volume Fraction of A/C'
    fAfC.scaler = Dictionary['fAfC']
    fAfC.datatype="Computational"
    
    PHASE = Property()
    PHASE.name = 'Phase'
    PHASE.scaler = Dictionary['PHASE']
    PHASE.datatype="Computational"

    
    REFERENCE = Reference(Dictionary['DOI'])
    
    PERSON = Person(Dictionary['NAME'],Dictionary['EMAIL'])
    
    chemical_system.classifications = Dictionary['CLASSIFICATION']
    chemical_system.properties = [SMILES,PHASE, RESOLUTION,FreeEnergy,CS,SPACEGroup,B,FreeEnergy,CHI,CELL,fAfC,NSC]
    chemical_system.references = REFERENCE
    chemical_system.person = PERSON
    return pif.dumps(chemical_system, indent=4)



def LAZY_maker(infile,outfile,DOI,SOFTWARE,METHOD,NAME,EMAIL):
    op = open(infile,'r')
    text = op.read().splitlines()
    op.close()
    DataDictionary = OutputParser(text)
    DataDictionary['DOI'] = DOI
    DataDictionary['METHOD'] = METHOD
    DataDictionary['NAME'] = NAME
    DataDictionary['EMAIL'] = EMAIL
    DataDictionary['fAfC'] = (DataDictionary['ARCH'][0][3]/np.sum(DataDictionary['ARCH'][0][3]))[0]
    if len(DataDictionary['ARCH'])==1:
        DataDictionary['NSC'] = 0
    else:
        DataDictionary['NSC'] = DataDictionary['ARCH'][1][3][0]
    
    DataDictionary['PHASE'] =infile.split('/')[-2][:-5]
    chainlist = list_convert_object(DataDictionary['ARCH'])
    G = molecule(chainlist).create_network()
    smiles = write_smiles(G)   
    DataDictionary['SMILES'] = smiles
    DataDictionary['CLASSIFICATION'] = ['Polymers','Bottlebrushes','Computational','SCFT']
    pif_dict = dict_to_pif(DataDictionary)
    print(f'Writing...{outfile}')
    op = open(f'{outfile}.json','w+')
    op.write(pif_dict)
    op.close()
#####################################################################################################
#Example of how to use
#####################################################################################################

# os.chdir('/home/tquah/Projects/TESTPHASES')
# #there are duplicates of 3D only true one is labeled by phase
# listoffiles = ['AGYRPhase/AGYR.out','ACYLPhase/ACYL.out','DISPhase/DIS.out','LAMPhase/LAM.out']


# DOI = '10.1021/acsmacrolett.0c00380'
# SOFTWARE = 'PolyFTS'
# METHOD = 'SCFT'
# NAME = 'Timothy Quah'
# EMAIL = 'timothy_quah@ucsb.edu'



# count = 0
# for file in listoffiles:
#     start = time.time()
#     outfile = f'test{count}'
#     LAZY_maker(file,outfile,DOI,SOFTWARE,METHOD,NAME,EMAIL)
#     count+=1
    
    
    
    
    
    

    


    
    # break
    # if i==3:
        

    # * A simulation cell has been initialized in 3 dimensions
    






