#!/usr/bin/env python3

import numpy as np
import pdb

'''
    Script to initialize fields when trying to generate SCFT seeds
    The idea is to use the idea that rho_A \approx w_- so that a reasonable initial guess
    for the w fields can be achieved using the density

    ...we'll see if this works
'''

# parameters
dim = 2
npw = [32,32]
gaussian_width = 0.2
centersA = [[0,0],[0.5,0.5]]
centersC = [[0,0.5],[0.5,0]]
magnitude = 10

# sanity checks
assert (len(npw) == dim)
for c in centersA:
    assert (len(c) == dim)
for c in centersC:
    assert (len(c) == dim)

# preliminaries
npw = tuple(npw)
inv_width2 = 1.0 / gaussian_width / gaussian_width
x = np.linspace(0,1,npw[0])
y = np.linspace(0,1,npw[1])

# w_plus \approx rhoA, w_1 \approx rhoC
w_A = np.zeros(npw)
w_C = np.zeros(npw)
w_B = np.random.random(npw)

# compute w_plus
for c in centersA:
  for ix,myx in enumerate(x):
    dx = (myx - c[0])
    #apply pbc
    if dx > 0.5:   dx -= 1.0
    if dx < -0.5:  dx += 1.0

    for iy,myy in enumerate(y):
      dy = (myy - c[1])
      if dy > 0.5:   dy -= 1.0
      if dy < -0.5:  dy += 1.0
      dr2 = dx*dx + dy*dy
      
      w_A[ix,iy] += np.exp(-0.5*dr2*inv_width2)
# set magnitude  [0,1] * magnitude
w_A = (w_A - np.min(w_A)) / (np.max(w_A) - np.min(w_A)) * magnitude
# set average to zero
w_A -= np.average(w_A)

# compute w_1 (rhoC)
for c in centersC:
  for ix,myx in enumerate(x):
    dx = (myx - c[0])
    #apply pbc
    if dx > 0.5:   dx -= 1.0
    if dx < -0.5:  dx += 1.0

    for iy,myy in enumerate(y):
      dy = (myy - c[1])
      if dy > 0.5:   dy -= 1.0
      if dy < -0.5:  dy += 1.0
      dr2 = dx*dx + dy*dy
      
      w_C[ix,iy] += np.exp(-0.5*dr2*inv_width2)

# set magnitude  [0,1] * magnitude
w_C = (w_C - np.min(w_C)) / (np.max(w_C) - np.min(w_C)) * magnitude
# set average to zero
w_C -= np.average(w_C)


# plot
import matplotlib.pyplot as plt 
xx, yy = np.meshgrid(x,y)
mesh = plt.pcolormesh(xx,yy,w_1)
plt.colorbar(mesh)
plt.savefig('fig_w_1.pdf')
plt.clf()
mesh = plt.pcolormesh(xx,yy,w_plus)
plt.colorbar(mesh)
plt.savefig('fig_w_plus.pdf')

# write file
filename = 'fields_test.in'
with open(filename,'w') as f:
  f.write("# Format version 3\n")
  f.write("# nfields = 3\n")
  f.write("# NDim = 2\n")
  f.write("# PW grid = 32 32 \n")
  f.write("# k-space data = 0 , complex data = 0\n")
  f.write("# Columns: x y fielddata\n")
  
  for ix,myx in enumerate(x):
    for iy,myy in enumerate(y):
      f.write(f"{ix} {iy} {w_A[ix][iy]} {w_B[ix][iy]} {w_C[ix][iy]}\n")
    f.write("\n")
    




