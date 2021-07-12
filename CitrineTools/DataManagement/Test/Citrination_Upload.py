#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 09:04:52 2021

@author: tquah
"""

from citrination_client import CitrinationClient
import os
import glob
apikey_path = os.path.join('/home/tquah/.citrine_api_key')
op = open(apikey_path,'r')
apikey = op.read().split('\n')
op.close()

client = CitrinationClient(apikey[0])


data_client = client.data


# file_path = "./test.json"
# dataset_id = 195077
# data_client.upload(dataset_id, file_path)

os.chdir('/home/tquah/Projects/TESTPHASES')

filelist = glob.glob('test*')

for i in range(len(filelist)):
    dataset_id = 51213104
    data_client.upload(dataset_id, filelist[i])
