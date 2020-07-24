#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 23 11:04:55 2020

@author: tquah
"""
import os
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


path = "/home/tquah/IMPORT_BRAID/"

files = os.listdir(path)
vtk_files = []
for file in files:
    if file.find('.vtk')>=0:
        vtk_files.append(os.path.join(path,file))        


filename = vtk_files[3]
import numpy
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN

reader = vtkStructuredPointsReader()
reader.SetFileName(filename)
reader.ReadAllVectorsOn()
reader.ReadAllScalarsOn()
reader.Update()

data = reader.GetOutput().GetPoints().GetData()
# nodes_nummpy_array = vtk_to_numpy(data)
