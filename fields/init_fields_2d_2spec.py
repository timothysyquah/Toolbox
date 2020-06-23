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
centers = [[0,0]]
magnitude = 10

# sanity checks
assert (len(npw) == dim)
for c in centers:
    assert (len(c) == dim)

# preliminaries
npw = tuple(npw)
inv_width2 = 1.0 / gaussian_width / gaussian_width
x = np.linspace(0,1,npw[0])
y = np.linspace(0,1,npw[1])

w_minus = np.zeros(npw)

# compute w_minus
for c in centers:
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
      
      w_minus[ix,iy] = np.exp(-0.5*dr2*inv_width2)

# set magnitude  [0,1] * magnitude
w_minus = (w_minus - np.min(w_minus)) / (np.max(w_minus) - np.min(w_minus)) * magnitude

# set average to zero
w_minus -= np.average(w_minus)

w_plus = np.random.random(npw)

# plot
import matplotlib.pyplot as plt 
xx, yy = np.meshgrid(x,y)
mesh = plt.pcolormesh(xx,yy,w_minus)
plt.colorbar(mesh)
plt.savefig('fig_w_minus.pdf')

# write file
filename = 'fields_test.in'
with open(filename,'w') as f:
  f.write("# Format version 3\n")
  f.write("# nfields = 2\n")
  f.write("# NDim = 2\n")
  f.write("# PW grid = 32 32 \n")
  f.write("# k-space data = 0 , complex data = 0\n")
  f.write("# Columns: x y fielddata\n")
  
  for ix,myx in enumerate(x):
    for iy,myy in enumerate(y):
      f.write(f"{ix} {iy} {w_plus[ix][iy]} {w_minus[ix][iy]}\n")
    f.write("\n")
    




