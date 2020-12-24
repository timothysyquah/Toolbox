#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  4 23:09:28 2020

@author: tquah
"""

import os
import execnet

def call_python_version(Version, Module, Function, ArgumentList):
    gw      = execnet.makegateway("popen//python=python%s" % Version)
    channel = gw.remote_exec("""
        from %s import %s as the_function
        channel.send(the_function(*channel.receive()))
    """ % (Module, Function))
    channel.send(ArgumentList)
    return channel.receive()


IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
toolpath = op.read()
op.close()

fieldcheckpath =toolpath+'/PolyFTSIO'

Fields_1 = '/home/tquah/IMPORT_BRAID/fieldstest/HEX/fields.in'
Fields_2 = '/home/tquah/IMPORT_BRAID/fieldstest/HEX/fields_k.bin'

def field_checker(fields_1,fields_2,tol):
    os.chdir(fieldcheckpath)
    result = call_python_version("2.7", "Fields_Check", "fields_compare",  
                                 [Fields_1, Fields_2]) 
    for point in result:
        if point<tol:
            print('Warning:Desired Field May Not Have been Formed')
            print(result)
    os.chdir(IDIR)


field_checker(Fields_1,Fields_2,0.99)