#!/usr/bin/env python3

import numpy as np
import glob
import os
import re
import pdb

def ExtractGapsFromOutput(filename):
  gaps = []
  with open(filename,'r') as f:
    for line in f:
      if 'Gap' in line: 
        l=line.split(' ')
        gap = float(re.sub('%','',l[-1]))
        gaps.append(gap)
  return gaps

def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    import re
    return [ atoi(c) for c in re.split('(\d+)', text) ]

import argparse as ap
parser = ap.ArgumentParser(description='This script plots the effect of varying the dielectric contrast between the A and C blocks')
parser.add_argument('--epsB',default=1.0,type=float,help='value of epsilon B to grab')
parser.add_argument('--inputdir',default='.',type=str,help='Input directory to where epsA*/epsB*/epsC* stack is located')
parser.add_argument('--output',type=str,help='name of output figure to generate')
parser.add_argument('--nosmooth',action='store_true',help='turn off smoothing of curves')
# Parse the command-line arguments
args=parser.parse_args()



#dir_to_analyze="delta0.00/fAfC0.25/ADIAPhase/"
#mpb_dir = dir_to_analyze + "/mpb/"
mpb_dir = args.inputdir
epsB=args.epsB
outfilename=args.output
idir=os.getcwd()

phase=""
for i,d in enumerate(mpb_dir.split('/')):
  if 'Phase' in d:
    if 'ADIA' in d:
      phase = "Diamond"
    elif 'AGYR' in d:
      phase = "Gyroid"
    
    
    
os.chdir(mpb_dir)

dirs = glob.glob('epsA*/epsB%0.1f/epsC*/' % epsB)
dirs.sort(key=natural_keys)
#print(dirs)
ndirs=len(dirs)
assert(ndirs != 0), "Error, no epsA*/epsB%0.1f/epsC* directories found. Try a different --inputdir or --epsB arguement" % epsB


#xx = np.zeros((N,N))
#yy = np.zeros((N,N))
epsAlist = []
epsClist = []
data = np.full(ndirs,0)

for i,d in np.ndenumerate(dirs):
  epsA = float(re.sub('[a-z,A-Z]','',d.split('/')[0]))
  epsB = float(re.sub('[a-z,A-Z]','',d.split('/')[1]))
  epsC = float(re.sub('[a-z,A-Z]','',d.split('/')[2]))
  if epsA not in epsAlist:
    epsAlist.append(epsA)
  if epsC not in epsClist:
    epsClist.append(epsC)

  mpb_outfile=f"{d}/mpb.out"
  if not os.path.isfile(mpb_outfile):
    print(f"{mpb_outfile} does not exist. Skipping")
    continue
  gaps = ExtractGapsFromOutput(mpb_outfile)
 
  #assert(len(gaps) <= 1)
  #if len(gaps) == 1:
  if len(gaps) > 1:
    print(f"Warning: multiple gaps found in {d}: {gaps}. Taking maximum.")
    data[i] = max(gaps)
  elif len(gaps) == 1:
    data[i] = gaps[0] 

# change directory back to root 
os.chdir(idir)

epsAlist = np.sqrt(epsAlist)
epsClist = np.sqrt(epsClist)

nepsA = len(epsAlist)
nepsC = len(epsClist)
depsA = epsAlist[1] - epsAlist[0]
depsC = epsClist[1] - epsClist[0]
assert (nepsA*nepsC == ndirs)

data = np.reshape(data,(nepsA,nepsC))

xx,yy = np.meshgrid(epsAlist,epsClist)
xx /= np.sqrt(epsB)
yy /= np.sqrt(epsB)
xmin,xmax=np.min(xx),np.max(xx)
ymin,ymax=np.min(yy),np.max(yy)


ndir_sqrt=np.sqrt(ndirs)
assert (ndir_sqrt == np.round(ndir_sqrt)) 
ndir_sqrt = np.int(ndir_sqrt)

# now plot
import matplotlib.pyplot as plt
import matplotlib as mpl

fig = plt.figure()
ax = fig.add_subplot(111)


pts = np.array([[1,1],np.sqrt([13,13])])
ax.plot(pts[:,0],pts[:,1],ls='--',color='k')
shift=1.5
#ax.text(6.5-shift,6.5+shift,'Double %s ($\epsilon_A = \epsilon_C$)' % phase ,va='center',rotation=45)
# plot single gyroid line
pts = np.array([[1,1],np.sqrt([13,1])])
ax.plot(pts[:,0],pts[:,1],ls='--',color='k')
shift=0.2
#ax.text(6.5,1-shift,'Single %s ($\epsilon_C = \epsilon_B$)' % phase,ha='center',va='top',rotation=0)
pts = np.array([[1,1],np.sqrt([1,13])])
ax.plot(pts[:,0],pts[:,1],ls='--',color='k')
#ax.text(1-shift,6.5,'Single %s ($\epsilon_A = \epsilon_B$)' % phase,va='center',ha='right',rotation=90)

# default colormap
#im=ax.imshow(data,cmap=mpl.cm.viridis,interpolation='None',origin='lower',extent=(min(epsAlist)-0.5,max(epsAlist)+0.5,min(epsClist)-0.5,max(epsClist)+0.5))

cbar_min=0
cbar_max=24
cbar_spacing = 2
ncbar_levels = int((cbar_max - cbar_min) / cbar_spacing) #10
data[data > cbar_max] = cbar_max
# fancy colormap
#colors = [(31/256.,120/256.,180/256.0,i) for i in np.linspace(0,1,ncbar_levels)] # blue
#colors = [(126/256.,30/256.,156/256.0,i) for i in np.linspace(0,1,ncbar_levels)] # purple
# even fancier colormap
colormap = plt.get_cmap('viridis') # using viridis
colors = []
for i,k in enumerate(np.linspace(0.1,1.00,ncbar_levels)):
  if i == 0: alpha = 0.0 #custom transparency
  else:
    alpha = 3*k # 1.0
    if alpha > 1: alpha = 1
  colors.append((colormap(k)[:3] + tuple([alpha])))
#colors=[(colormap(k)[:3] + tuple([k])) for k in np.linspace(0, 1.00, ncbar_levels)] # override alpha in 4th column
cmap = mpl.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=ncbar_levels)

if args.nosmooth:
  #xx -= 0.5*depsA
  #yy -= 0.5*depsC
  im=ax.pcolormesh(xx,yy,data,cmap=cmap,vmin=cbar_min,vmax=cbar_max,shading='flat')
else:
  #N=100
  #x = np.array(epsAlist)
  #y = np.array(epsClist)
  #xnew = np.linspace(xmin,xmax,N)
  #ynew = np.linspace(ymin,ymax,N)
  #xxnew, yynew = np.meshgrid(xnew,ynew)
  #from scipy import interpolate
  #f = interpolate.interp2d(x,y,data,kind='cubic')
  #datanew = f(xnew,ynew)
  #im=ax.pcolormesh(xxnew,yynew,datanew,cmap=cmap,vmin=cbar_min,vmax=cbar_max,shading='flat')

  im=ax.contourf(xx,yy,data,levels=ncbar_levels,cmap=cmap,vmin=cbar_min,vmax=cbar_max)
  
  # original interpolation
  #im=ax.imshow(data,interpolation='gaussian',cmap=cmap,origin='lower',extent=(xmin,xmax,ymin,ymax),vmin=cbar_min,vmax=cbar_max)


cbar= fig.colorbar(im)
cbar.set_ticks(np.linspace(cbar_min,cbar_max,ncbar_levels+1))

#ax.set_xlabel('$\epsilon_A / \epsilon_B$')
#ax.set_ylabel('$\epsilon_C / \epsilon_B$')
ax.set_aspect('equal','box')
ax.set_xlim([0.9,xmax])
ax.set_ylim([0.9,ymax])
ax.set_xlabel('Index contrast, $n_C / n_B$')
ax.set_ylabel('Index contrast, $n_A / n_B$')
cbar.set_label('Complete-gap size (%)')


if outfilename:
  plt.savefig(outfilename)
plt.show()

