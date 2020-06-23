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
ndim = 3
npw = [32,32,32]
nspecies = 3

phase='C15'
if phase == 'C15':
  gaussian_width = 0.1
  # set centersA, centerC using C15 positions
  general_positions = [[0,0,0],[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]] # C15, space group 227 (Fd-3m)  
  special_positions_8a = [[0,0,0],[0.75,0.25,0.75]] # Wyckoff position 8a
  special_positions_16d = [[5./8,5./8,5./8],[3./8,7./8,1./8],[7./8,1./8,3./8],[1./8,3./8,7./8]]
  centersA = []
  centersC = []
  for gp in general_positions:
    for sp in special_positions_8a:
      pos = np.add(gp,sp)
      for i in range(ndim):
        while (pos[i] >= 1): pos[i] -= 1
        while (pos[i] <  0): pos[i] += 1
      centersA.append(list(pos))

    for sp in special_positions_16d:
      pos = np.add(gp,sp)
      for i in range(ndim):
        while (pos[i] >= 1): pos[i] -= 1
        while (pos[i] <  0): pos[i] += 1
      centersC.append(list(pos))


#centersC = [[0,0.5],[0.5,0]]
magnitude = 10

# sanity checks
assert (len(npw) == ndim)
for c in centersA:
    assert (len(c) == ndim)
for c in centersC:
    assert (len(c) == ndim)

# preliminaries
npw = tuple(npw)
inv_width2 = 1.0 / gaussian_width / gaussian_width
x = np.linspace(0,1,npw[0])
y = np.linspace(0,1,npw[1])
z = np.linspace(0,1,npw[2])
xx,yy,zz = np.meshgrid(x,y,z)


def update_species_field_gaussians(field, centers):
  # compute w_plus
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

        for iz,myz in enumerate(z):
          dz = (myz - c[2])
          if dz > 0.5:   dz -= 1.0
          if dz < -0.5:  dz += 1.0

          dr2 = dx*dx + dy*dy + dz*dz

          field[ix,iy,iz] -= np.exp(-0.5*dr2*inv_width2)
  # set magnitude  [0,1] * magnitude
  field = (field - np.min(field)) / (np.max(field) - np.min(field)) * magnitude
  # set average to zero
  field -= np.average(field)

def update_species_field_levelset(phase,field,xx,yy,zz,originshift=(0,0,0)):
  ''' 
  Reference: Wohlgemuth, M., Yufa, N., Hoffman, J. & Thomas, E. L. Macromolecules 34, 6083–6089 (2001).
  '''
  pi = np.pi

  # Fig 1a (Eq 5)
  #F = np.cos(2*np.pi*x) + np.cos(2*np.pi*y) + np.cos(2*np.pi*z)
  # Fig 2A gyroid
  #F =  np.sin(2*pi*y) * np.cos(2*pi*z) \
  #      +  np.sin(2*pi*z) * np.cos(2*pi*x) \
  #      +  np.sin(2*pi*x) * np.cos(2*pi*y) 

  # Fig 2d, double gyroid 
  #F = 0.8*(np.sin(4*pi*x)*np.sin(2*pi*z)*np.cos(2*pi*y) \
  #         + np.sin(4*pi*y)*np.sin(2*pi*x)*np.cos(2*pi*z) \
  #         + np.sin(4*pi*z)*np.sin(2*pi*y)*np.cos(2*pi*x))\
  #   - 0.2*(np.cos(4*pi*x)*np.cos(4*pi*y) \
  #          + np.cos(4*pi*y)*np.cos(4*pi*z) \
  #          + np.cos(4*pi*z)*np.cos(4*pi*x))

  if phase == 'diamond':
    # Fig 3a, diamond (Eq 13)
    x = xx + originshift[0] #0.125
    y = yy + originshift[1] # 0.125
    z = zz + originshift[2] # 0.125
    field -= np.cos(2*pi*x)*np.cos(2*pi*y)*np.cos(2*pi*z) \
        + np.sin(2*pi*x)*np.sin(2*pi*y)*np.cos(2*pi*z)\
        + np.sin(2*pi*x)*np.cos(2*pi*y)*np.sin(2*pi*z)\
        + np.cos(2*pi*x)*np.sin(2*pi*y)*np.sin(2*pi*z)
  else:
    raise RuntimeError(f"Invalid phase {phase}")

  # Fig 3c, double diamond
  #x += 0.25
  #y += 0.25
  #z += 0.25
  #F = 0.5*(np.sin(2*pi*x)*np.sin(2*pi*y) \
  #         + np.sin(2*pi*y)*np.sin(2*pi*z) \
  #         + np.sin(2*pi*x)*np.sin(2*pi*z)) \
  #    + 0.5*np.cos(2*pi*x)*np.cos(2*pi*y)*np.cos(2*pi*z) 
  #return F 

  # set magnitude  [0,1] * magnitude
  field = (field - np.min(field)) / (np.max(field) - np.min(field)) * magnitude
  # set average to zero
  field -= np.average(field)



# w_plus \approx rhoA, w_1 \approx rhoC
w_A = np.zeros(npw)
w_C = np.zeros(npw)
w_B = np.random.random(npw)

if phase == 'C15':
  # set fields using centers
  update_species_field_gaussians(w_A,centersA)
  update_species_field_gaussians(w_C,centersC)
elif phase == 'diamond':
  # set diamond fields, note that the different origin shifts are what give the alternating network
  update_species_field_levelset("diamond",w_A,xx,yy,zz,originshift=(0.125,0.125,0.125))
  update_species_field_levelset("diamond",w_C,xx,yy,zz,originshift=(0.625,0.625,0.625))


# A is matrix to map from species fields to exchange fields

# I got this from PolyFTS in ModelMolecularAuxFieldBase.cpp in variable m_adLinearMap_SpeciesFieldsToHSFields

# this was for Model_BlockPolymerMelt, chiN12=19, chiN13 = 72, chiN23 = 19 
# will need to regenerate for different chi values
Ainv = np.matrix(np.zeros((nspecies,nspecies)))
Ainv[0,:] = [0.89241734306928799, 0.45121091053868922, -1.3436282536079771]
Ainv[1,:] = [-0.45121091053868922, 0.89241734306928799, -0.44120643253059877]
Ainv[2,:] = [0, 0, 1]

# this was for Model_BlockPolymerMelt, chiN12=26, chiN13 = 70, chiN23 = 26 
#Ainv[0,:] = [0.87526325903077462, 0.48364680024872186, -1.3589100592794965]
#Ainv[1,:] = [-0.48364680024872186, 0.87526325903077462, -0.39161645878205276]
#Ainv[2,:] = [0, 0, 1]



w_exch = np.zeros((nspecies,npw[0],npw[1],npw[2]))
for ix,myx in enumerate(x):
  for iy,myy in enumerate(y):
    for iz,myz in enumerate(z):
      w_spec_tmp = np.matrix([w_A[ix,iy,iz], w_B[ix,iy,iz], w_C[ix,iy,iz]]).T
      w_exch_tmp = np.matmul(Ainv, w_spec_tmp)
      for ispec in range(nspecies):
        w_exch[ispec,ix,iy,iz] = w_exch_tmp[ispec]




# write file
def write_file(filename,fieldtype):
  with open(filename,'w') as f:
    f.write("# Format version 3\n")
    f.write(f"# nfields = {nspecies}\n")
    f.write("# NDim = 3\n")
    f.write(f"# PW grid = {npw[0]} {npw[1]} {npw[2]} \n")
    f.write("# k-space data = 0 , complex data = 0\n")
    f.write("# Columns: x y z fielddata\n")
    
    for iz,myz in enumerate(z):
      for iy,myy in enumerate(y):
        for ix,myx in enumerate(x):
          if fieldtype == 'species':
            f.write(f"{ix} {iy} {iz} {w_A[ix][iy][iz]} {w_B[ix][iy][iz]} {w_C[ix][iy][iz]}\n")
          elif fieldtype == 'exchange':
            f.write(f"{ix} {iy} {iz} {w_exch[0][ix][iy][iz]} {w_exch[1][ix][iy][iz]} {w_exch[2][ix][iy][iz]}\n")
          else:
            raise RuntimeError(f"invalid fieldtype {fieldtype}")
      f.write("\n")
    

write_file('species_fields_generated.in','species')
write_file('exchange_fields_generated.in','exchange')


