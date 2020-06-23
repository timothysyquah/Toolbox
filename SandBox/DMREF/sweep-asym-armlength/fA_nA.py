#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 28 11:10:47 2020

@author: tquah
"""

import numpy as np

nTot = 100
NA = 10
NB = 40

def fA_calc(nA,nB,NA,NB):
    return (nA*(NA+1))/((nA*(NA+1))+(nB*(NB+1)))
def nA_calc(fA,nTot,NA,NB):
    return (fA*nTot*(NB+1))/((1-fA)*(NA+1)+fA*(NB+1))



fA = np.arange(0.01,0.5+1e-6,0.04)
nA_check = nA_calc(fA,nTot,NA,NB)
nA_int = np.array(list(set(list(np.around(nA_check))+list([nTot]))))
del_loc = np.where(nA_int==nTot)[0]
nA_int = np.delete(nA_int,del_loc)
nA_int = np.sort(nA_int)
fA_act = fA_calc(nA_int,nTot-nA_int,NA,NB)

print(fA_act)
print(nA_int)