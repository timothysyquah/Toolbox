#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 12:27:59 2021

@author: tquah
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from stats import *


# plt.close('all')
path = '/home/tquah/Projects/TESTAVERAGE/compare_runtimes/operators10000.dat'








# path = '/home/tquah/Projects/TESTAVERAGE/compare_runtimes/operators10000.dat'
# df = pd.read_csv(path,sep = " ")
# header = list(df)


# desiredheader = header[1:-1]
# df = df.drop(columns = [header[-1]])
# df = df.drop(columns = [header[-2]])
# df.columns = desiredheader

# parameter = 'StressXX.Real'


# df['100'] = df[parameter].rolling(100, min_periods=1).mean()
# df['10'] = df['StressXX.Real'].rolling(10, min_periods=1).mean()
# df['1000'] = df['StressXX.Real'].rolling(1000, min_periods=1).mean()
# df['CMA'] = df['StressXX.Real'].expanding().mean()
# # df['EMA_0.1'] = df['StressXX.Real'].ewm(alpha=0.1, adjust=False).mean()
# # df['EMA_0.3'] = df['StressXX.Real'].ewm(alpha=0.05, adjust=False).mean()


# plt.figure()
# # plt.plot(df['step'],df['StressXX.Real'])
# # plt.plot(df['step'],df['10'])
# # plt.plot(df['step'],df['100'])
# plt.plot(df['step'],df['1000'])
# # plt.plot(df['step'],df['CMA'])
# # plt.plot(df['step'],df['EMA_0.1'])
# # plt.plot(df['step'],df['EMA_0.3'])

# path = '/home/tquah/Projects/TESTAVERAGE/compare_runtimes/operators5000.dat'
# df = pd.read_csv(path,sep = " ")
# header = list(df)


# desiredheader = header[1:-1]
# df = df.drop(columns = [header[-1]])
# df = df.drop(columns = [header[-2]])
# df.columns = desiredheader

# parameter = 'StressXX.Real'


# df['100'] = df[parameter].rolling(100, min_periods=1).mean()
# df['10'] = df['StressXX.Real'].rolling(10, min_periods=1).mean()
# df['1000'] = df['StressXX.Real'].rolling(1000, min_periods=1).mean()
# df['CMA'] = df['StressXX.Real'].expanding().mean()
# # df['EMA_0.1'] = df['StressXX.Real'].ewm(alpha=0.1, adjust=False).mean()
# # df['EMA_0.3'] = df['StressXX.Real'].ewm(alpha=0.05, adjust=False).mean()


# # plt.plot(df['step'],df['StressXX.Real'])
# # plt.plot(df['step'],df['10'])
# # plt.plot(df['step'],df['100'])
# plt.plot(df['step'],df['1000'])
# # plt.plot(df['step'],df['CMA'])
# # plt.plot(df['step'],df['EMA_0.1'])
# # plt.plot(df['step'],df['EMA_0.3'])