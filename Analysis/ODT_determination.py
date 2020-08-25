#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 12:14:32 2020

@author: tquah
"""

import numpy as np
Nbb = 100
Nsc = np.array([1,5,15,20])
alpha = Nsc/Nbb
chiNeff = 24.9*alpha**(0.839)+10.5