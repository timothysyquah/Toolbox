#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 24 13:23:46 2020

@author: tquah
"""
import glob
import argparse



nodes = initialize_nodes(args.dirs, args.filename,args.keywrd)
boundaryholder = calc_phase_boundaries(nodes)
boundaryholder.plot(args.outfig,args.plottype, nodes=nodes,xlabel=args.xlabel, ylabel=args.ylabel,axisrange=args.axisrange,n=args.dim,aspect=args.aspect,refPhase=args.refphase)
if args.raw != '':
    print("Saving free energy curve data to \'%s\'" % args.raw)
    boundaryholder.write(args.raw,dim=1)
