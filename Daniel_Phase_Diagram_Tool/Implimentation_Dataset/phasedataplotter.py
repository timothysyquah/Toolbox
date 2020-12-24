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

def findPhaseBoundary(fAidx,epsidx,phaseidx):
	assert(len(phaseidx) == 2)
	result = []
	for fPBidx0 in (0,3):
		for fPBidx1 in (1,2):
			tmpBdy = FindEdgeIntersection([[fA[fAidx+shifts[fPBidx0][0]], eps[epsidx+shifts[fPBidx0][1]]], [fA[fAidx+shifts[fPBidx1][0]], eps[epsidx+shifts[fPBidx1][1]]]],
																	phaseidx)
			if len(tmpBdy) > 0:
				result.append(tmpBdy)

	return result

def FindEdgeIntersection(pts,phases):
	# given two points (for which we have energies) on an edge, and two phases, find where their free energies intersect
	assert(len(phases) == 2)
	assert(len(pts) == 2)

	dE = []
	for fidx in range(2): # loop over points
		idxlist = []
		for pidx in range(2):
			tmpidx = np.where(np.logical_and(data[phases[pidx]][:,0]==pts[fidx][0],data[phases[pidx]][:,1]==pts[fidx][1] ))[0]
			if (tmpidx.size == 0):
				print('Error: Missing point for phaseidx ' + phase_names[phases[m]] + ' at fA=' + str(pts[fidx][0]) + ', eps=' + str(pts[fidx][1]))
				return []
			else:
				idxlist.append(tmpidx[0])
		dE.append(data[phases[0]][idxlist[0],2] - data[phases[1]][idxlist[1],2])


	result = []
	if (np.abs(dE[0]) < 1e-8  and np.abs(dE[1]) < 1e-8 ):
		result = [] # if the energy diff is zero for both points, there is no phase boundary
	elif (np.abs(dE[0]) < 1e-8 ):
		result = pts[0] 
	elif (np.abs(dE[1]) < 1e-8 ):
		result = pts[1] 
	elif (dE[0]*dE[1] < 0):
		t = dE[0]/(dE[0] -dE[1])
		result = [t*pts[1][0] + (1-t)*pts[0][0], t*pts[1][1] + (1-t)*pts[0][1]]
	else:
		result = []
	return result

def InterpolatePhaseEnergy(pts,phase,ipt):
	if not ipt:
		print("Error: interpolation point is empty")
		return []
	# given two points (for which we have energies), interpolate the energy of a phase at a third colinear point
	# probably need to check for co-linearity
	t = np.linalg.norm(np.array(ipt)-np.array(pts[0]))/np.linalg.norm(np.array(pts[1])-np.array(pts[0]))
	idxlist= []
	for fidx in range(2):
		tmpidx = np.where(np.logical_and(data[phase][:,0]==pts[fidx][0],data[phase][:,1]==pts[fidx][1] ))[0]
		if (tmpidx.size == 0):
			print('Error: Missing point for phaseidx ' + phase_names[phases[m]] + ' at fA=' + str(pts[fidx][0]) + ', eps=' + str(pts[fidx][1]))
			return np.inf
		else:
			idxlist.append(tmpidx[0])
	result = (1-t)*data[phase][idxlist[0],2] + t* data[phase][idxlist[1],2]
	return result
	
	

	

			
draw_boundaries = True
marker = ('+', 'o', '*','v','^','<','>','s','p','h','H','x')
#markers=[',', '+', '-', '.', 'o', '*']
ymin = -0.1
ymax = 0.05
fmin = 0.01
fmax = 0.5

phase_names = ['BCC','DIS','HEX','A15','SIGMA']
nphases = len(phase_names)
data = [] # store data as a list of np arrays

for j in range(len(phase_names)):
	data.append(np.genfromtxt(phase_names[j] + '.dat',skip_header=1))
# 	data[-1][:,0] /= 50.
# 	data[-1][np.lexsort((data[-1][:,0],data[-1][:,1]))]
	#print(np.lexsort((data[-1][:,0],data[-1][:,1])))
	#print(data[-1])


fig, ax = plt.subplots(1,sharex='all')
fig.set_size_inches(10,5)

# now compute minimum energy phase
# first need a list of all the fA samples
fA_set = set(data[0][:,0])
for j in range(1,len(data)):
	fA_set.union(set(data[j][:,0]))
fA = np.sort(np.fromiter(fA_set,float,len(fA_set)))
# all values of fA now in the set
# also need a list of all the eps samples
eps_set = set(data[0][:,1])
for j in range(1,len(data)):
	eps_set.union(set(data[j][:,1]))
eps = np.sort(np.fromiter(eps_set,float,len(eps_set)))


#dfA = 0.5
#deps = 0.25
#fA_interp, eps_interp = np.mgrid[min(fA):max(fA):dfA],np.mgrid[min(eps):max(eps):deps]
#tmp1, tmp2 = np.mgrid[min(fA):max(fA):dfA,min(eps):max(eps):deps]
##fA_interp, eps_interp = np.mgrid([np.minimum(fA):np.maximum(fA):dfA,min(eps):max(eps):deps])
#data_interp = []
#for j in range(len(data)):
#	data_interp.append(sp.interpolate.griddata((data[j][:,0],data[j][:,1]),data[j][:,2],(tmp1,tmp2),method='linear'))
#
#convex_hull = []
#for j in range(len(fA_interp)):
#	for k in range(len(eps_interp)):
#		Emin = 0.
#		idxmin = phase_names.index('DIS')
#		Emin = data_interp[idxmin][j,k]
#		for l in range(len(data_interp)):
#			if (data_interp[l][j,k] < Emin ):
#				idxmin = l
#				Emin = data_interp[l][j,k]
#		convex_hull.append([fA_interp[j],eps_interp[k],idxmin,Emin])
#		#Emin_tot = np.minimum(Emin, Emin_tot)
#print(convex_hull)
# all values of eps now in the set
# now find the minimum free energy at each value of fA
#Emin_tot = 0.
convex_hull = []
for j in range(len(fA)):
	for k in range(len(eps)):
		Emin = 0.
		idxmin = phase_names.index('DIS')
		tmp = np.where(np.logical_and(data[idxmin][:,0]==fA[j],data[idxmin][:,1]==eps[k]))[0]
		if (tmp.shape[0] != 0):
			Emin = data[idxmin][tmp[0],2]
		else:
			Emin = np.inf
		for l in range(len(data)):
			tmp = np.where(np.logical_and(data[l][:,0]==fA[j],data[l][:,1]==eps[k]))[0]
			if ( tmp.shape[0] != 0 ):
				idx = tmp[0]
				if ( Emin - data[l][idx,2] > 1e-8 ):
					idxmin = l
					Emin = data[l][idx,2]
		convex_hull.append([fA[j],eps[k],idxmin,Emin])
		#Emin_tot = np.minimum(Emin, Emin_tot)
#print(convex_hull)
scatter_data = [];
for j in range(len(phase_names)):
	scatter_data.append(np.zeros((0,2)))
for j in range(len(convex_hull)):
	idx = convex_hull[j][2]
	scatter_data[idx] = np.append(scatter_data[idx],[[convex_hull[j][0],convex_hull[j][1]]],axis=0)
nstablephases = 0 # keep track of how many phases are stable
stable_phase_names = []
for j in range(nphases):
	if (scatter_data[j].size != 0 ):
		nstablephases += 1
		stable_phase_names.append(phase_names[j])
		ax.scatter(scatter_data[j][:,0],scatter_data[j][:,1])

ax.legend(stable_phase_names,bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
#plt.show()

# Now find the actual phase boundaries via linear interpolation between nearest neighbors
# On a regular 2D grid there are Ny*(Nx-1) + Nx*(Ny -1) nearest neighbor pairs there are also
# Nphases * (Nphases -1 ) /2 possible phase pairings for a possible boundary
#nphasepairs = int(nstablephases * (nstablephases - 1) /2)
#
## strategy for finding points on phase boundaries:
## 1. start in bottom left corner of phase diagram (smallest fA and eps)
## 2. Define rectangle using the bottom leftmost 4 point
## 3. Determine how many different phases are stable on the corners of the rectangle
##		a. If 1 -> do nothing. There is only 1 stable phase over this region and no phase boundaries
##		b. If 2 -> Determine which edges of the square connect nodes with different stable phases
##			 Do Linear interpolation along these edges to find the phase boundary points
##		c. If 3 or more -> Do a superset of 2. May Have to truncate some lines

if draw_boundaries:
	phase_boundaries = [];
	shifts=[[0,0],[1,0],[0,1],[1,1]]
	for j in range(len(fA) - 1):
		for k in range(len(eps) -1):
			phases=[]
			for l in range(len(shifts)):
				tmp = [ x for x in convex_hull if (x[0] == fA[j + shifts[l][0]] and  x[1] == eps[k + shifts[l][1]])] 
				phases.append( tmp[0][2] )
			phases_norepeats = list(set(phases))
			#print(fA[j],eps[k],phases,)
			if len(phases_norepeats) > 1:
				bdypts = []
				# determine convex hull along each edge
				for l in (0,3):
					for m in (1,2):
						# does this edge connect two different phases?
						#print(phases,l,m,phases[l],phases[m])
						midpts = [];
						if ( phases[l] != phases[m] ):
							endpts = [[fA[j+shifts[l][0]],eps[k + shifts[l][1]]],[fA[j + shifts[m][0]],eps[k + shifts[m][1]]]]
							xm = FindEdgeIntersection(endpts,[phases[l],phases[m]])
							if len(xm) == 0:
								print('Error: Empty point on first edge pass')
							midpts.append([xm,InterpolatePhaseEnergy(endpts,phases[l],xm),[phases[l],phases[m]]])
							for n in range(len(data)): # loop over all other phases
								delidx_list = []
								add_list = []
								if (not  np.any(n == np.array([phases[l],phases[m]]))):
									for o in range(len(midpts)): # loop over all intermediate points
										if (midpts[o][1] - InterpolatePhaseEnergy(endpts,n,midpts[o][0]) > 1e-8 ): # This phase is lower energy, we should consider it
											delidx_list.append(o) # this point no longer matter
											# create two new intersection points since this phase goes in the middle
											for p in range(2): # loop over the two new points we are adding
												xmtmp = FindEdgeIntersection(endpts,[midpts[o][2][p],n])
												if len(xmtmp) == 0:
													print('Error: No intersection found when trying to construct edge boundaries')
													print('Phases: ' + phase_names[n] + ' ' + phase_names[midpts[o][2][p]] )
													print('Edge Endpoints: ')
													print(endpts)
													print('End point phases: ' + phase_names[phases[l]] + ' ' + phase_names[phases[m]])
												add_list.append([xmtmp,InterpolatePhaseEnergy(endpts,n,xmtmp),[midpts[o][2][p],n]])
								# delete the necessary old points
								for o in sorted(delidx_list, reverse=True):
									del midpts[o]
								# add in the new points
								midpts += add_list
							
							#print(endpts,midpts)
							bdypts += midpts
				#print(bdypts)
				#WARNING: This next section is specialized assuming three or fewer phases exist in this tile
				# for an even number of boundaries, then each boundary point should have a mate
				# add these pairs to the boundary list, then remove them
				if (len(bdypts)%2 == 0):
					while (len(bdypts) > 0):
						for l in range(1,len(bdypts)):
							if (sorted(bdypts[0][2]) == sorted(bdypts[l][2])):
								phase_boundaries.append([bdypts[0][0],bdypts[l][0]])
								del bdypts[l]
								del bdypts[0]
								break
				else:
					bdytmp = [];
					slopevec =[]
					# compute all three boundaries, then find intersectio
					for l in range(len(bdypts)):
						bdytmp.append(findPhaseBoundary(j,k,bdypts[l][2]))
						#print(bdytmp[l])
						slopevec.append((bdytmp[l][1][1] - bdytmp[l][0][1])/(bdytmp[l][1][0] - bdytmp[l][0][0]))
					
					# strictly all three boundaries should intersect at one point, but they'll probably be a little off. Instead use average of all three individual intersections
					favg = 0
					eavg = 0
					for l in range(len(bdypts)-1):
						for m in range(l+1,len(bdypts)):
							fcross=-((bdytmp[l][0][1] - slopevec[l]*bdytmp[l][0][0]) - (bdytmp[m][0][1] - slopevec[m]*bdytmp[m][0][0]))/(slopevec[l]-slopevec[m])
							ecross=bdytmp[l][0][1] + slopevec[l]*(fcross - bdytmp[l][0][0] )
							favg += fcross
							eavg += ecross
					favg /= (len(bdypts)*(len(bdypts)-1))/2.
					eavg /= (len(bdypts)*(len(bdypts)-1))/2.
					#print(favg, eavg)
					#print(bdypts)
					for l in range(len(bdypts)):
						phase_boundaries.append([bdypts[l][0],[favg,eavg]])	
	
				# Now we actually compute boundaries based on phases that have a boundary along the edge
				#for newidx in range(len(bdypts)):
								
	
	for j in range(len(phase_boundaries)):
		ax.plot([phase_boundaries[j][0][0],phase_boundaries[j][1][0]],[phase_boundaries[j][0][1],phase_boundaries[j][1][1]],c='k')

ax.set_xlabel('volume fraction of A, $\phi_A$')
ax.set_ylabel('conformational asymmetry, $\epsilon$')


plt.savefig('Phase_diagram.pdf')
plt.show()


