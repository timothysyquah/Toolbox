#!/usr/bin/env python3

import os
import execnet
import sys
# sys.path.append('/home/tquah/toolbox_github/lib/iotools')
# from iotools import ReadBinFile


def call_python_version(Version, Module, Function, ArgumentList):
    gw      = execnet.makegateway("popen//python=python%s" % Version)
    channel = gw.remote_exec("""
        from %s import %s as the_function
        channel.send(the_function(*channel.receive()))
    """ % (Module, Function))
    channel.send(ArgumentList)
    return channel.receive()


path_to_tools = '/home/tquah/.timtools'
op = open(path_to_tools,'r')
toolpath = op.read()
op.close()

fieldcheckpath =toolpath+'/PolyFTSIO'


def field_checker(fields_1,fields_2):
    os.chdir(fieldcheckpath)
    result = call_python_version("2.7", "Fields_Check", "fields_compare",  
                                 [fields_1, fields_2]) 
    return result


field_2 = '/home/tquah/toolbox_github/readdensity/GYR_Fields_2/density.bin'
field_1 = '/home/tquah/toolbox_github/readdensity/GYR_Fields/density.bin'
result = field_checker(field_1,field_2)

