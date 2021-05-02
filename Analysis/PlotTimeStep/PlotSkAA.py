#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 15:16:30 2021

@author: tquah
"""
import os
import numpy as np
import matplotlib.pyplot as plt
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_Test/rotation_test/Nbb40/SKAA.dat'
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/prototypedir/nbb109.0/LAMPhase_CL_noncubic_test/SKAA.dat'
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/replicant_test/c_8_4/SKAA.dat'
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DomainSpacing/c_9.0/SKAA.dat'
# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/c_8.0/SKAA.dat'
# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/c_6.2/SKAA.dat'
IDIR = os.getcwd()
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_6.2'
os.chdir(path)
data = np.loadtxt('SKAA.dat')
plt.close('all')
fig = plt.figure()
ax = plt.gca()
# ax.plot(data[:,0],data[:,1],'or',label = 'L=62 l')
# ax.plot(data[:,0],data[:,1],'r')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
# os.chdir(IDIR)


# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_7.0'
# os.chdir(path)
# data = np.loadtxt('SKAA.dat')
# ax.plot(data[:,0],data[:,1],'og',label = 'L=70 l')
# ax.plot(data[:,0],data[:,1],'g')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
# os.chdir(IDIR)

# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_8.0'
# os.chdir(path)
# data = np.loadtxt('SKAA.dat')
# ax.plot(data[:,0],data[:,1],'ob',label = 'L=80 l')
# ax.plot(data[:,0],data[:,1],'b')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
# os.chdir(IDIR)
# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_9.2'
# os.chdir(path)
# data = np.loadtxt('SKAA.dat')
# ax.plot(data[:,0],data[:,1],'ok',label = 'L=92 l')
# ax.plot(data[:,0],data[:,1],'k')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
# os.chdir(IDIR)



path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/extrarefine/c_9.2_1000'
os.chdir(path)
data = np.loadtxt('SKAA.dat')
ax.plot(data[:,0],data[:,1],'^r',label = 'L=92 l')
ax.plot(data[:,0],data[:,1],'r')
plt.xlim(0,1.0)
plt.ylim(-1,200)
ax.set_xlabel(r'$q(b/\sqrt{6})$')
ax.set_ylabel(r'$S_{AA}(q)/CN$')
ax.set_yscale('symlog')
plt.tight_layout()
os.chdir(IDIR)


# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/extrarefine/c_9.2_2000'
# os.chdir(path)
# data = np.loadtxt('SKAA.dat')
# ax.plot(data[:,0],data[:,1],'^b',label = 'L=92 l')
# ax.plot(data[:,0],data[:,1],'b')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
# os.chdir(IDIR)


# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/extrarefine/c_9.2_3000'
# os.chdir(path)
# data = np.loadtxt('SKAA.dat')
# ax.plot(data[:,0],data[:,1],'^m',label = 'L=92 l')
# ax.plot(data[:,0],data[:,1],'m')
# plt.xlim(0,1.5)
# plt.ylim(-1,200)
# ax.set_xlabel(r'$q(b/\sqrt{6})$')
# ax.set_ylabel(r'$S_{AA}(q)/CN$')
# ax.set_yscale('symlog')
# plt.tight_layout()
os.chdir(IDIR)


ax.legend()

plt.savefig('/home/tquah/Figures/skaa1.pdf',dpi = 300)
