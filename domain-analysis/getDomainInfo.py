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

import pdb
import scipy.spatial #import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mpl

np.set_printoptions(linewidth=240)
# assume 3D for now

infile="density.bin"
coords, fields = io.ReadBinFile(infile)
# need to get box sizes to account for PBC, or just get delta and NPW from coords

domainanalyzer = DomainAnalyzer(coords,fields)
domainanalyzer.setDensityThreshold(0.5)
ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=False,add_periodic_domains=False)

with open("domains.dat","w") as outfile:
	outfile.write("x y z area volume IQ\n")
	for j in range(ndomains):
		outfile.write('{0: e} {1: e} {2: e} {3} {4} {5}\n'.format(com[j,0],com[j,1],com[j,2],area[j],vol[j],IQ[j]))
#print(np.hstack((com,area[:,np.newaxis],vol[:,np.newaxis],IQ[:,np.newaxis])))
#print(ndomains)
#print(com)
#print(area)
#print(vol)
#print(IQ)

