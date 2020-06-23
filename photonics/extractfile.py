#!/usr/bin/env python3

import os 
import shutil
keyword ='density'
filetype = '.bin'
# This is to get the directory that the program  
# is currently running in. 
dir_path = os.path.dirname(os.path.realpath(__file__)) 
newdirectory_path = dir_path+'/'+keyword+'_export'
if os.path.isdir(newdirectory_path):    
    os.rmdir(dir_path+'/'+keyword+'_export')
os.mkdir(dir_path+'/'+keyword+'_export')

uscore = '_'
for root, dirs, files in os.walk(dir_path): 
    for file in files:  
        if root!=newdirectory_path:
            # change the extension from '.mp3' to  
            # the one of your choice. 
            if file.endswith(keyword+filetype): 
                print(root)
                input_path = root+'/'+str(file)
                component = input_path.split('/')
                output_name = component[-1]+uscore+\
                                component[-4]+uscore+\
                                component[-3]+uscore+\
                                component[-2]
                newpath = os.path.join(newdirectory_path,output_name)
                shutil.copyfile(input_path,newpath)
