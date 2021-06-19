#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 10:30:49 2021

@author: tquah

Here I try to reproduce the method found in "Fluctuations, Phase Transitions, and Latent Heat in Short Diblock
Copolymers: Comparison of Experiment, Simulation, and Theory" by Timothy M. Gillard, Pavani Medapuram, David C. Morse,* and Frank S. Bates*

https://pubs.acs.org/doi/pdf/10.1021/acs.macromol.5b00277

"""
import numpy as np


phiA = np.array([0.51,0.53,0.51,0.43,0.63,0.52,0.62,0.61,0.69])
b_A = 9.2
b_B = 8.2
v0 = 118

cn= ((b_A*phiA+b_B*(1-phiA))**3)/v0/6/np.sqrt(6)
