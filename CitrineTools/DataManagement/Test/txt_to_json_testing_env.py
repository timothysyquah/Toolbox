#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 20:11:34 2021

@author: tquah
"""

import glob
from txt_to_json import * 
import json
import re

    


path = '/home/tquah/Projects/CITRINE/**/*.in'
filelist = glob.glob(path)



final_filelist = []
for i in range(len(filelist)):
    if 'field' not in filelist[i]:
        final_filelist.append(filelist[i])
        

# first step


for file in final_filelist:
    op = open(file,'r+')
    txtfile = op.read()

    myString = re.sub(r"[\n\t]*", "", txtfile)

    op.close()