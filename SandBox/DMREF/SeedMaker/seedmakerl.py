#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 17:42:31 2020

@author: tquah
"""

from make_input_file import make_input
import numpy as np
import os
import pdb
import shutil
import time
import argparse
from itertools import combinations
import subprocess
from runner_general import *



if __name__ == '__main__':


    IDIR = os.getcwd()

    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-PolyFTS', '--polyFTS_path', action='store', default='/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release',help='PolyFts Path',type = str)
    # parser.add_argument('-s', '--sweeptype', action='store', default='chi',help='Sweep Flory-Huggins (chi), armlength (Nsc), or backbone (Nbb) or Volume Fraction (f)',type = str)
    parser.add_argument('-stat', '--chain_stat', action='store', default='DGC',help='Chain Statistics', type = str)
    parser.add_argument('-p', '--phase', action='store',nargs='+', default=['DIS','LAM'],help='Phases that should be investigated',type = str)
    # parser.add_argument('-sg', '--spacegroup', action='store',nargs='+', default=[None,None],help='Spacegroups for Phases',type = spacegroup_lgic)
    parser.add_argument('-nt', '--NThreads', action='store',nargs='+', default=[1,1],help='Number of Threads',type = int)
    parser.add_argument('-sp', '--seed_path', action='store', default=os.path.join(IDIR,'SEEDS'),help='Path to Seeds',type = str)
    #Simulation Details
    parser.add_argument('-ds', '--contour_step', action='store', default=0.1,help='Contour step for CGC',type = float)
    parser.add_argument('-cl', '--chain_label', action='store', default='AB-Bottlebrush',help='Chain Label',type = str)
    parser.add_argument('-dt', '--timestep', action='store', default=0.01,help='Timestep size',type = float)
    parser.add_argument('-nref', '--reference_length', action='store', default=1,help='Reference Length',type = int)
    parser.add_argument('-ss', '--stressscale', action='store', default=0.1,help='Stress Mobility',type = float)
    parser.add_argument('-fs', '--forcescale', action='store',nargs='+', default=[1.0,1.0],help='Force Mobility',type = float)
    parser.add_argument('-L0', '--initial_box_size', action='store',nargs='+', default=[5.0,5.0],help='Initial Box Size',type = float)
    parser.add_argument('-npw', '--number_of_planewaves', action='store',nargs='+', default=[512,512],help='Number of Planewaves',type = int)
    parser.add_argument('-cu', '--cell_updater', action='store', default='Broyden',help='Cell updater usually Broyden or Euler',type = str)
    parser.add_argument('-fu', '--field_updater', action='store', default='EMPEC',help='Field Updater EMPEC or SIS',type = str)
    parser.add_argument('-ftol', '--force_tol', action='store', default=1e-5,help='Force stopping tolerance',type = float)
    parser.add_argument('-stol', '--stress_tol', action='store', default=1e-4,help='stress stopping tolerance',type = float)
    
    #Chain Archetecture
    parser.add_argument('-bsp', '--backbone_species', action='store',nargs='+', default=[1,2],help='Backbone Species',type=int)
    parser.add_argument('-sas', '--sidearm_species', action='store',nargs='+', default=[[1],[2]],help='Sidearm Species',type=sidearm_species_lgic)
    parser.add_argument('-sac', '--sidearm_composition', action='store',nargs='+', default=[[1.0],[1.0]],help='Sidearm Composition',type=sidearm_composition_lgic)
    parser.add_argument('-space', '--sidearm_spacing', action='store',nargs='+', default=[1,1],help='Sidearm Spacing',type=int)


    parser.add_argument('-nbb', '--backbone', action='store', default=99.0,help='Length of Polymer Backbone',type=float)



    #sweepable parameters
    parser.add_argument('-nbbmin', '--backbone_min', action='store', default=99.0,help='Length of Polymer Backbone min',type=float)
    parser.add_argument('-nbbmax', '--backbone_max', action='store', default=99.0,help='Length of Polymer Backbone max',type=float)
    parser.add_argument('-dnbb', '--backbone_step', action='store', default=1.0,help='Length of Polymer Backbone step',type=float)
    parser.add_argument('-nsc', '--sidechains', action='store',nargs='+', default=[20,20],help='Length of Side Chains ',type = float)


    parser.add_argument('-nscmin', '--sidechains_min', action='store',nargs='+', default=[20,20],help='Length of Side Chains min',type = float)
    
    parser.add_argument('-nscmax', '--sidechains_max', action='store',nargs='+', default=[20,20],help='Length of Side Chains max',type = float)
    parser.add_argument('-ndsc', '--sidechains_step', action='store',nargs='+', default=[20,20],help='Length of Side Chains step',type = float)

    parser.add_argument('-chimin', '--floryhuggins_interaction_parameter_min', action='store', default=[0.1],help='Flory-Huggins Interaction Parameter min (12 13 23)',type = float)
    parser.add_argument('-chimax', '--floryhuggins_interaction_parameter_max', action='store', default=[0.1],help='Flory-Huggins Interaction Parameter max (12 13 23)',type = float)
    parser.add_argument('-dchi', '--floryhuggins_interaction_parameter_step', action='store', default=[-0.005],help='Flory-Huggins Interaction Parameter step (12 13 23)',type = float)

    parser.add_argument('-fmin', '--volumefraction_min', action='store',nargs='+', default=[0.5,0.5],help='Volume Fraction min',type = float)
    parser.add_argument('-fmax', '--volumefraction_max', action='store',nargs='+', default=[0.5,0.5],help='Volume Fraction max',type = float)
    parser.add_argument('-df', '--volumefraction_step', action='store',nargs='+', default=[-0.02,0.02],help='Volume Fraction step',type = float)
    parser.add_argument('-dirs', '--directory_structure', action='store',nargs='+', default=['chi','nbb','nsc','f'],help='Order directory for sweepable parameters',type = str)
    parser.add_argument('-dirnl', '--directory_nameloc', action='store',nargs='+', default=[0,None,0,0],help='Location for directory name location',type = directory_nameloc_lgic)

