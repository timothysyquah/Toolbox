#!/usr/bin/env python3
import sys
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
import itertools
np.set_printoptions(linewidth=180)
plt.rcParams.update({
	'font.size': 16, 
	'figure.autolayout': True, 
	'lines.linewidth': 2, 
	'axes.linewidth':2, 
	'xtick.major.width': 2,
	'xtick.major.size': 8,
	'xtick.minor.width': 1,
	'xtick.minor.size': 4,
	'ytick.major.width': 2,
	'ytick.major.size': 8,
	'ytick.minor.width': 1,
	'ytick.minor.size': 4,
	'legend.handletextpad': 0.05,
	'legend.columnspacing': 0.05,
	'lines.markersize': 10
	})
#markers=[',', '+', '-', '.', 'o', '*']
xmin = 0.1
xmax = 0.35
phase_key = ()

cololist = ['m','k','y','g','b','b','r']
#colodict = {
#		"HEX": [0.8,0.7,1.0],
#		"DIS": [0.7,0.7,0.7],
#		"BCC": [1,0.95,0.65],
#		"A15": [0.6,0.9,0.6],
#		"FCC": [0.6,0.7,1],
#		"HCP": [0.6,0.7,1],
#		"SIGMA": [1,0.6,0.7]
#}
colodict = {
    "HEX": [151,118,181],
    "DIS": [178,178,178],
    "BCC": [255,190,64],
    "A15": [102,163,133],
    "FCC": [67,124,155],
    "HCP": [67,124,155],
    "SIGMA": [255,85,119],
    "GYR": [124,155,67],
    "LAM": [0,255,0],
}


phasewidth = []
phase_order = []
phase_key = []
sysnames = ['Standard','SCdiff','BBdiff','Graftdiff']
namestoplot = ['    All $b=1$','$b_{B,sc}=3/2$', '$b_{B,bb}=2/3$','$z_A=1/6$']
epstoplot = ['$\epsilon \\approx 2$', '$\epsilon \\approx 2$ ','$\epsilon \\approx 3$', '$\epsilon \\approx 3$ ']


nsys = len(sysnames)
for sysidx in range(len(sysnames)):
	phase_names = ()
	data = [] # store data as a list of np arrays
	
	with open(sysnames[sysidx] + '.dat','r') as f:
		next(f) # first line is header
		for line in f:
			if (line[0] == "#" ): # add the name of the phase
				phase_names += (line[1:].rstrip(),)
				data.append(np.empty( shape=(0,2) ))
			elif (line != "\n"): # add the data for that phase
				tmp = np.fromstring(line.rstrip(),sep=' ')
				data[-1] = np.vstack((data[-1],tmp[np.newaxis,:]))
	if ( sysidx == 0):
		for j in range(len(data)):
			data[j][:,0] /= 600
	if ( sysidx == 1 or sysidx == 2):
		for j in range(len(data)):
			data[j][:,0] /= 200
	if ( sysidx == 3):
		for j in range(len(data)):
			data[j][:,0] *= 8./6./1200.

	#for j in range(len(data)): # loop over phases
	#	for k in range(data[j].shape[0]):
	#		idx=np.where(fAmap[:,0] == data[j][k,0])[0][0]
	#		if ( sysidx < 3 ):
	#			data[j][k,0] = fAmap[idx,2]
	#		else:
	#			data[j][k,0] = fAmap[idx,4]
	#for j in range(len(data)):
	#	data[j][:,0] /= 600
	# now compute minimum energy phase
	# first need a list of all the fA samples
	fA_set = set(data[0][:,0])
	for j in range(1,len(data)):
		fA_set = fA_set.union(set(data[j][:,0]))
	fA = np.sort(np.fromiter(fA_set,float,len(fA_set)))
	# all values of fA now in the set
	# now find the minimum free energy at each value of fA
	Emin_tot = 0.
	convex_hull = []
	for j in range(len(fA)):
		idxmin = phase_names.index('DIS')
		idx = np.where(data[idxmin][:,0]==fA[j])[0][0]
		Emin = data[idxmin][idx,1]
		for l in range(len(data)):
			tmp = np.where(data[l][:,0]==fA[j])[0]
			if (tmp.size != 0):
				idx = tmp[0]
				if ((Emin - data[l][idx,1]) > 1e-6  ):
					idxmin = l
					Emin = data[l][idx,1]
		convex_hull.append((fA[j],idxmin,phase_names[idxmin],Emin))
		Emin_tot = np.minimum(Emin, Emin_tot)
	#print(convex_hull)

	#print(fA)
	
	# now find boundaries via linear interpolation
	blist = []
	for j in range(1,len(convex_hull)):
		if (convex_hull[j][1] != convex_hull[j-1][1]):
			(p1, p2) = (convex_hull[j-1][1], convex_hull[j][1])
			(f0, f1) = (convex_hull[j-1][0], convex_hull[j][0])
			dE0 = data[p1][np.where(data[p1][:,0]==f0)[0][0],1] - data[p2][np.where(data[p2][:,0]==f0)[0][0],1]
			dE1 = data[p1][np.where(data[p1][:,0]==f1)[0][0],1] - data[p2][np.where(data[p2][:,0]==f1)[0][0],1]
			f = dE0*(f1-f0)/(dE0-dE1) + f0
			blist.append((f,convex_hull[j-1][2],convex_hull[j][2]))
	
	#print(blist)
	
	# create ordered list of phase and stability range thickness
	if (len(phasewidth) ==0):
		firstpass = True
		for j in range(len(phase_names)):
			phasewidth.append(np.zeros(shape=(nsys)))
			phase_key.append(phase_names[j])
	else:
		firstpass = False
	
	idx = phase_key.index(blist[0][1])
	if firstpass:
		phase_order.append(idx)
	phasewidth[idx][sysidx] = blist[0][0] - xmin
	for j in range(1,len(blist)):
		#print(np.where(phase_names==blist[j][1]))
		idx = phase_key.index(blist[j][1])
		if firstpass:
			phase_order.append(idx)
		phasewidth[idx][sysidx] = blist[j][0] - blist[j-1][0]
	idx = phase_key.index(blist[j][2])
	if firstpass:
		phase_order.append(idx)
	phasewidth[idx][sysidx] = np.minimum(xmax,fA.max()) - blist[-1][0]



#print(phasewidth)
print(phase_order)
phase_order = [1,2,4,3,0]
print(phase_key)

#for j in range(len(phasewidth)):
# phasewidth[j][5] = 0.
# phasewidth[j][4] = 0.


fig, ax = plt.subplots()
ax.invert_yaxis()
ax.set_xlim([xmin,xmax])
ax.set_xlabel('$f_A$')


curleft = np.ones(nsys)*xmin
for j in range(0,len(phase_order)):
	ax.barh(namestoplot,phasewidth[phase_order[j]],align='center',left=curleft,label=phase_key[phase_order[j]],color=np.array(colodict[phase_key[phase_order[j]]])/255)
	#ax.text(curleft[0] + phasewidth[phase_order[j]][0]/2.,0,phase_key[phase_order[j]],ha='center',va='center')
	curleft += phasewidth[phase_order[j]]

ax2 = ax.twinx()
ax2.invert_yaxis()

ax2.barh(epstoplot,np.zeros(len(epstoplot)))

#ax.axhline(1.5,c='k',lw=4)

#ax.text(xmin+ phasewidth[phase_order[0]][0]/2., 0,'DIS',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][0]+phasewidth[phase_order[1]][0]/2., 0,'BCC',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][0]+phasewidth[phase_order[1]][0]+phasewidth[phase_order[2]][0]/2., 0,'$\sigma$',ha='center',va='center')
#ax.text(xmin+phasewidth[phase_order[0]][0]+phasewidth[phase_order[1]][0]+phasewidth[phase_order[2]][0] + phasewidth[phase_order[3]][0]/2., 0,'$\sigma$',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][0]+phasewidth[phase_order[1]][0]+phasewidth[phase_order[2]][0] + phasewidth[phase_order[3]][0]+ phasewidth[phase_order[4]][0]/2., 0,'HEX',ha='center',va='center')

ax.text(xmin+phasewidth[phase_order[0]][1]+phasewidth[phase_order[1]][1]/2., 1,'BCC',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][1]+phasewidth[phase_order[1]][1]+phasewidth[phase_order[2]][1]/2., 1,'$\sigma$',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][1]+phasewidth[phase_order[1]][1]+phasewidth[phase_order[2]][1] + phasewidth[phase_order[3]][1]+ phasewidth[phase_order[4]][1]/2., 1,'HEX',ha='center',va='center')
#
ax.text(xmin+phasewidth[phase_order[0]][2]+phasewidth[phase_order[1]][2]/2., 2,'BCC',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][2]+phasewidth[phase_order[1]][2]+phasewidth[phase_order[2]][2]/2., 2,'$\sigma$',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][2]+phasewidth[phase_order[1]][2]+phasewidth[phase_order[2]][2] + phasewidth[phase_order[3]][2]/2., 2,'A15',ha='center',va='center')
#
ax.text(xmin+phasewidth[phase_order[0]][3]+phasewidth[phase_order[1]][3]/2., 3,'BCC',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][3]+phasewidth[phase_order[1]][3]+phasewidth[phase_order[2]][3]/2., 3,'$\sigma$',ha='center',va='center')
ax.text(xmin+phasewidth[phase_order[0]][3]+phasewidth[phase_order[1]][3]+phasewidth[phase_order[2]][3] + phasewidth[phase_order[3]][3]/2., 3,'A15',ha='center',va='center')

#ax.legend(loc='lower right')
plt.savefig('CaseStudies.pdf')
plt.show()
