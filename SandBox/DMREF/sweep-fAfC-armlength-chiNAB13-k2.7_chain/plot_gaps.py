#!/usr/bin/env python3

import numpy as np
import glob
import os
import re

def ExtractGapsFromOutput(filename):
  gaps = []
  with open(filename,'r') as f:
    for line in f:
      if 'Gap' in line: 
        l=line.split(' ')
        gap = float(re.sub('%','',l[-1]))
        gaps.append(gap)
  return gaps

#plottype='Nsc'
plottype='fAfC'
if plottype == 'Nsc':
  dirs = glob.glob('Nsc*/fAfC0.25/')
elif plottype == 'fAfC':
  dirs = glob.glob('Nsc0.040/fAfC*/')
eps = [7.8,1,1]

phases = ['AGYR32','ADIA32', 'ACYL']
ncol = len(phases) + 1
ndir = len(dirs)
data = np.full((ndir,ncol),np.nan)
for i,d in enumerate(dirs):
  if plottype=='Nsc':
    chiNAB = float(re.sub('Nsc','',d.split('/')[0]))
    data[i][0] = chiNAB
  elif plottype=='fAfC':
    fAfC = float(re.sub('fAfC','',d.split('/')[1]))
    data[i][0] = fAfC

  for iphase,phase in enumerate(phases):
    mpb_outfile=f"{d}/{phase}Phase/mpb-{eps[0]:.1f}-{eps[1]:.1f}-{eps[2]:.1f}.out"
    if not os.path.isfile(mpb_outfile):
      print(f"{mpb_outfile} does not exist. Skipping")
      continue
    
    gaps = ExtractGapsFromOutput(mpb_outfile)
    #assert(len(gaps) <= 1)
    #if len(gaps) == 1:
    if len(gaps) > 1:
      print(f"Warning: multiple gaps found in {d}: {gaps}. Taking maximum.")
      data[i][iphase+1] = max(gaps)
    elif len(gaps) == 1:
      data[i][iphase+1] = gaps[0] 
    elif len(gaps) == 0:
      data[i][iphase+1] = 0.0
 
data = data[data[:,0].argsort()] # sort by 0th column
print(data) 
np.savetxt(f"gaps_{plottype}.dat",data, header=f"{plottype} {phases}")

import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
if plottype=='fAfC':
  ax.set_xlabel(r'$\mathcal{f}_{\rm A} = \mathcal{f}_{\rm C} $')
ax.set_ylabel('Gap size (%)')
for iphase,phase in enumerate(phases):
  ax.plot(data[:,0],data[:,iphase+1],label=f"{phase}",marker='.')
ax.legend()
plt.savefig(f'fig_gaps_{plottype}.pdf')
plt.show()

