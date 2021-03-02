#!/usr/bin/env python3

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)

import fieldtools 
import iotools as io
import numpy as np
from domaintools import DomainAnalyzer
import regex as re
import pdb
import scipy.spatial #import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

def make_grid(dirlist,rules,key_locations):
    xlist = []
    ylist = []
    for i in range(len(dirlist)):
        listofstr = re.findall("\d+\.\d+", dirlist[i])
        listofvalue = np.array([float(i) for i in listofstr])
        x = rules[0](listofvalue[key_locations[0]])
        y = rules[1](listofvalue[key_locations[1]])
        xlist.append(x)
        ylist.append(y)
    return np.vstack((xlist,ylist)).transpose()


phases = ['BCC','A15','SIGMA']
IDIR = os.getcwd()
key_locations = [3,[1,2]]
#rules
tuplerule = lambda x: np.sqrt((x[1]+1)/(x[0]+1))
fA = lambda x : x*1

rules = [fA,tuplerule]


for i in range(len(phases)):
    outfile = open(phases[i]+'_domain_analysis.dat','w')
    outfile.write("fA eps x y z area volume IQ\n")

    op = open(phases[i]+'.imp','r')
    listofdir = op.read().splitlines()
    array = make_grid(listofdir,rules,key_locations)
    op.close()
    os.chdir('/home/tquah/Projects/DMREF/sweep-asym-armlength_corrected')
    for j in range(len(listofdir)):
        infile = os.path.join(listofdir[j],'density.bin')
        coords, fields = io.ReadBinFile(infile)
        # need to get box sizes to account for PBC, or just get delta and NPW from coords
        
        domainanalyzer = DomainAnalyzer(coords,fields)
        domainanalyzer.setDensityThreshold(0.5)
        ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=False,add_periodic_domains=False)

        for k in range(ndomains):
            outfile.write('{0} {1} {2: e} {3: e} {4: e} {5} {6} {7}\n'.format(array[j,0], array[j,1],\
                                                                              com[k,0],com[k,1],com[k,2],\
                                                                                  area[k],vol[k],IQ[k]))
        outfile.close()
        # with open("domains.dat","w") as outfile:
        #  	outfile.write("x y z area volume IQ\n")
        #  	for j in range(ndomains):
        # 		outfile.write('{0: e} {1: e} {2: e} {3} {4} {5}\n'.format(com[j,0],com[j,1],com[j,2],area[j],vol[j],IQ[j]))

    
    
    os.chdir(IDIR)
    # for j in range(len())


# np.set_printoptions(linewidth=240)
# # assume 3D for now

# infile="density.bin"
# #print(np.hstack((com,area[:,np.newaxis],vol[:,np.newaxis],IQ[:,np.newaxis])))
# #print(ndomains)
# #print(com)
# #print(area)
# #print(vol)
# #print(IQ)

