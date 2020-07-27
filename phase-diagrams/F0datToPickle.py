#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 22:54:27 2020

@author: tquah
"""

import numpy as np
import argparse
import glob,re
import os,sys
import matplotlib.pyplot as plt
from phase_data_plotter import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-o', '--outfig', action='store', default='',help='name of output figure file')
    parser.add_argument('--raw', action='store', default='',help='name of raw output file (for plotting in another program')
    parser.add_argument('-t', '--plottype', action='store', default='simplecolors',help='type of plot to generate')
    parser.add_argument('--xlabel', action='store', default=r"$f_A$",help='label for xaxis, can use \'$\' to write latex')
    parser.add_argument('--zlabel', action='store', default=r"$\tau$")
    parser.add_argument('--ylabel', action='store', default=r"$\chi N$",help='')
    parser.add_argument('--axisrange', action='store', nargs=4, default=[None,None,None,None],help='')
    parser.add_argument('--linecutoff', action='store', nargs='+', default=[1e30,1e30],help='maximum length of lines to draw in phase diagrams, useful to clean them up')
    parser.add_argument('-n','--dim',action='store',default=None,help='Number of dimensions to plot phase data in \n   1 => Free energy curves\n   2 => Phase Diagram \n   3 => 3d phase diagram  (guesses by default)')
    parser.add_argument('-i','--interp_dimension',action='store',default=[0],nargs='+',help='Dimensions to interpolate the phase diagram along ex: [0,1] would interpolate in 2 dimensions')
    parser.add_argument('-p','--plotstyle3d',action='store',default='flat',help='This argument changes the 3d plot style. Flat => multiple graphs with different linestyles on top of each other')
    parser.add_argument('--stylesheet',action= 'store',default=os.path.dirname(os.path.realpath(sys.argv[0]))+'/better_style.mplstyle',help='This argument is the Matplotlib stylesheet that will be used for graphing') #the default is located in the directory this script is located at
    parser.add_argument('--aspect',action='store',default=None,help='The aspect ratio for the outputted figure use 1 for a square fig, works for a 2d graph right now')
    parser.add_argument('-r', '--refphase', action='store', default=None,help='name of phase to reference to, only matters if 1d')
    print("IMPLEMENT CUSTOM AXIS RANGES AND LABELS FROM COMMAND LINE")
    args = parser.parse_args()
    #fnmeIn="F0_phases.dat"
    #dirs=glob.glob("tau*/phiA*");
    if os.path.isfile(args.stylesheet):
        plt.style.use(args.stylesheet)
        print('Graphing using the {} stylesheet'.format(args.stylesheet))
    else:
        print('WARNING: No stylesheet found at {} graphing with default Matplotlib settings'.format(args.stylesheet))
        print(os.path.dirname(os.path.realpath(sys.argv[0])))

    # Guess how many dimensions you're plotting
    if not args.dim:
        #Here it finds the where there are numbers in the first directory string, and then checks if those numbers are differnt in subsequent directory strings
        dirs=args.dirs


        folders = re.split('/',dirs[0])
        print('Trying to guess how many dimensions to plot in')
        args.dim = 0
        locs = [] # stores the indicies of directories in tree that contain numbers
        nums1 = [] # stores the actual numbers in each subdirectory of tree
        names=[]
        index = 0
        for folder in folders:
           if '*' in folder:
                raise ValueError("There is a '*' in you directories. Looks like a wildcard didn't get expanded. Check your path, something is likely wrong!")
           if re.search('[0-9]',folder):
               locs.append(index)
               first_num = re.sub('[^0-9.]',"",re.split('_',folder)[-1])#get the last number in the folder, sometimes there are multiple numbers
               nums1.append(float(first_num))
           index +=1
        for i in range(len(nums1)):
           #in every path if the number at the current location is different add one to the dimension
           for mydir in dirs[1:]:
               if mydir[-1] == '/':
                 raise RuntimeError(f'Directory "{mydir}" contains a trailing slash, this messes up parsing. Please remove and try again.')

               mysubdir = re.split('/',mydir)[locs[i]]
               if float(re.sub('[^0-9.]',"",re.split('_',mysubdir)[-1])) != nums1[i]:#grab the last number following the underscore since there may be multiple numbers in a directory
                    args.dim += 1
                    names.append(re.sub('[0-9,.]','',mysubdir))
                    break
        if args.dim == 0 or args.dim >3:
            raise ValueError('could not guess how many dimensions display in please specify with the -n flag')

        print('Graphing in {0} dimensions according to guess. If this is incorrect specify the number manually with the -n flag'.format(str(args.dim)))

        if args.dim == 1:
            args.xlabel = names[0]
            args.ylabel = 'Free Energy'
        elif args.dim == 2:
            args.ylabel = names[0]
            args.xlabel = names[1]
    else:
        args.dim = int(args.dim)


    if args.axisrange != [None, None, None, None]:
        args.axisrange = [float(x) for x in args.axisrange]
    print('axis ranges:',args.axisrange)
    print('xlabel:',args.xlabel)
    print('ylabel:',args.ylabel)
    print('zlabel:',args.zlabel)

    if args.interp_dimension != [0]:
        args.interp_dimension = [int(i) for i in args.interp_dimension]


    if args.dim == 1:
        nodes = initialize_nodes(args.dirs, args.filename)
        boundaryholder = calc_phase_boundaries(nodes)
        boundaryholder.plot(args.outfig,args.plottype, nodes=nodes,xlabel=args.xlabel, ylabel=args.ylabel,axisrange=args.axisrange,n=args.dim,aspect=args.aspect,refPhase=args.refphase)
        if args.raw != '':
            print("Saving free energy curve data to \'%s\'" % args.raw)
            boundaryholder.write(args.raw,dim=1)
