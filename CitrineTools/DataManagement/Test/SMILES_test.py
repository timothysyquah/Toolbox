import numpy as np
import os
import networkx as nx
from pysmiles import write_smiles, fill_valence
import matplotlib.pyplot as plt
import matplotlib
# from rdkit import Chem
import scipy as sp

#matplotlib.use('TKAgg')
import string
plt.close('all')
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

#stolen from stackoverflow 
    


#backbone
# N_{bb} | Species Order | Composition
chainlist = []
# chainlist.append(chain('DGC',100,[1,2,3],[0.25,0.5,0.25]))
# chainlist.append(sidechain('DGC',5,[1],[1.0],1,25,25))
# chainlist.append(sidechain('DGC',5,[2],[1.0],26,75,50))
# chainlist.append(sidechain('DGC',5,[3],[1.0],76,100,25))


chainlist.append(chain('DGC',100,[1,2,3],[0.2,0.6,0.2]))
chainlist.append(sidechain('DGC',3,[1],[1.0],1,4,4))
chainlist.append(sidechain('DGC',3,[2],[1.0],5,16,12))
chainlist.append(sidechain('DGC',3,[3],[1.0],17,20,4))

G = molecule(chainlist).create_network()
nx.draw(G, with_labels=True, font_weight='bold')
nodelist = []
for i in range(1,len(G.nodes)+1,1):
    nodelist.append(G.nodes[i]['element'])
A = nx.adjacency_matrix(G).todense()

# Chem.MolToSmiles(MolFromGraphs(nodelist, A))

# N_sc | Species order | Composition|  Graft Start | Graft End | Number of sidechains

# sidechain_info_1 = [5,[1],[1.0],1,25,25]
# sidechain_info_2 = [5,[2],[1.0],26,75,50]
# sidechain_info_3 = [5,[3],[1.0],76,100,25]














#for idx, ele in enumerate(string):
#    mol.nodes[idx]['element'] = ele
#    mol.add_edges_from([(idx,idx+1)])
#    count+=1
    
    
    
print(write_smiles(G))
# [O-]C(=O)C([C])([C])[C]
# fill_valence(mol, respect_hcount=True)
#print(write_smiles(mol))
# [O-]C(=O)C(C)(C)C



#plt.subplot(121)

#nx.draw(G, with_labels=True, font_weight='bold')
#plt.subplot(122)

# nx.draw(G)
#plt.show()