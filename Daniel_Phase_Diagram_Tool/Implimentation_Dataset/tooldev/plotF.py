#!/usr/bin/env python3
import sys
import numpy as np
import scipy as sp
import scipy.stats
import matplotlib.pyplot as plt
import itertools
from scipy.spatial import Delaunay,Voronoi
import os
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

class Edge:
    def __init__(self,ptidx):
        # ptidx should be a list of two integer indices for the pts in pt_array to use
        phases = [None]*2
        self.endpts = [pt_array[ptidx[0],:], pt_array[ptidx[1],:]]
        for j in range(2):
            phases[j] = [ x for x in convex_hull if (x[0] == self.endpts[j][0] and  x[1] == self.endpts[j][1])][0][2]

        midpts = [];
        self.bdypts = []
        if ( phases[0] != phases[1]):
            xm = FindEdgeIntersection(ptidx,phases)
            #print(xm)
            if len(xm) == 0: # if xm returns empty, just assume the crossing is in the middle
                xm = [(self.endpts[0][0]+self.endpts[1][0])/2., (self.endpts[0][1]+self.endpts[1][1])/2.]
                E = (FindPhaseEnergy(self.endpts[0],phases[0]) + FindPhaseEnergy(self.endpts[1],phases[1]))/2.
                midpts.append([xm,E,[phases[0],phases[1]]])
            else:
                midpts.append([xm,InterpolatePhaseEnergy(self.endpts,phases[0],xm),[phases[0],phases[1]]])
            #print(midpts)
            
            for n in range(len(data)): # loop over all other phases
                delidx_list = []
                add_list = []
                if (not  np.any(n == np.array(phases))):
                    for o in range(len(midpts)): # loop over all intermediate points
                        if (midpts[o][1] - InterpolatePhaseEnergy(self.endpts,n,midpts[o][0]) > tol ): # This phase is lower energy, we should consider it
                            delidx_list.append(o) # this point no longer matter
                            # create two new intersection points since this phase goes in the middle
                            for p in range(2): # loop over the two new points we are adding
                                xmtmp = FindEdgeIntersection(ptidx,[midpts[o][2][p],n])# if there is no intersection, then phase n must be lower in energy than phase midpts[o][2][p] over the entire interval and we don't add anything to the list of addpts
                                # there is actually a potential problem here. The phases could also intersect at some point where a previously considered phase is lower energy
                                if len(xmtmp) > 0: 
                                    add_list.append([xmtmp,InterpolatePhaseEnergy(self.endpts,n,xmtmp),[midpts[o][2][p],n]])
                # delete the necessary old points
                for o in sorted(delidx_list, reverse=True):
                    del midpts[o]
                # add in the new points
                midpts += add_list
            
            #print(endpts,midpts)
            self.bdypts += midpts
            #print(self.bdypts)

    def getBdyPts(self):
        return self.bdypts


def InterpolatePhaseEnergy(pts,phase,ipt):
    #if (not ipt):
    #    print("Error: interpolation point is empty")
    #    return np.nan
    # given two points (for which we have energies), interpolate the energy of a phase at a third colinear point
    # probably need to check for co-linearity
    t = np.linalg.norm(np.array(ipt)-np.array(pts[0]))/np.linalg.norm(np.array(pts[1])-np.array(pts[0]))
    idxlist= []
    Elist = [None]*2
    for fidx in range(2):
        Elist[fidx] = FindPhaseEnergy(pts[fidx],phase)
    result = (1-t)*Elist[0] + t*Elist[1]
    return result
        
def FindPhaseEnergy(pt,phase):
    tmpidx = np.where(np.logical_and(data[phase][:,0]==pt[0],data[phase][:,1]==pt[1] ))[0]
    if (tmpidx.size == 0):
        #print('Missing point for phaseidx ' + phase_names[phase] + ' at fA=' + str(pt[0]) + ', eps=' + str(pt[1]))
        return np.nan
    else:
        return data[phase][tmpidx[0],2]

def FindEdgeIntersection(ptidx,phases):
    # given two points (for which we have energies) on an edge, and two phases, find where their free energies intersect
    assert(len(phases) == 2)
    assert(len(ptidx) == 2)

    dE = [None]*2
    for ipt in range(2): # loop over points
        Elist = []
        for jpt in range(2): # loop over phases
            Etmp = FindPhaseEnergy(pt_array[ptidx[ipt],:],phases[jpt])
            if np.isnan(Etmp): # if we are missing an energy, don't return anything
                return []
            else:
                Elist.append(Etmp)
        dE[ipt] = Elist[0] - Elist[1]

    result = []
    if (np.abs(dE[0]) < tol  and np.abs(dE[1]) < tol ):
        result = [] # if the energy diff is zero for both points, there is no phase boundary
    elif (np.abs(dE[0]) < tol ): #if the energy diff is below the tol at one or the other point, thats where the boundary is
        result = pt_array[ptidx[0],:]  
    elif (np.abs(dE[1]) < tol ):
        result = pt_array[ptidx[1],:]  
    elif (dE[0]*dE[1] < 0): #there must be a sign change in dE for there to be a boundary so that t is positive
        t = dE[0]/(dE[0] -dE[1])
        result = [(1-t)*pt_array[ptidx[0],0]  + t*pt_array[ptidx[1],0], (1-t)*pt_array[ptidx[0],1] + t*pt_array[ptidx[1],1] ]
    else:
        result = []
    return result
    

            
draw_boundaries = True
marker = ('+', 'o', '*','v','^','<','>','s','p','h','H','x')
color  = ('tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan')

tol = 1e-10
#markers=[',', '+', '-', '.', 'o', '*']
#ymin = -0.1
#ymax = 0.05
fmin = 0.05
fmax = 0.40

phase_names = ['BCC','DIS','HEX','A15','SIGMA']
# phase_names = ['BCC','DIS','HEX','LAM','GYR']

exportname = '/home/tquah/Figures/comparison_2.png'
mainarray = np.loadtxt('../mainarray.txt')

# os.chdir("/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected/PHASE_FREE_ENERGY")
# os.chdir("/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected/PHASE_FREE_ENERGY")
# keywords = ['over']
keywords = ['all','under']
# keywords = ['under']
keywords = ['all']

fig, ax = plt.subplots(1,sharex='all')
loc = np.where(mainarray[:,0]<0.4)[0]


counter_keyword = 0
for keyword in keywords:
    os.chdir(f"/home/tquah/toolbox/Daniel_Phase_Diagram_Tool/Implimentation_Dataset/PHASE_FREE_ENERGY_{keyword}")
    # os.chdir("/media/tquah/TOSHIBA EXT/Projects/sweep-asym-armlength_BCC_fix/PHASE_FREE_ENERGY")
    # phase_names = ['BCC','DIS','HEX','LAM','GYR']
    nphases = len(phase_names)
    data = [] # store data as a list of np arrays
    
    for j in range(len(phase_names)):
        data.append(np.genfromtxt(phase_names[j] + '.dat',skip_header=1))
        data[-1][np.lexsort((data[-1][:,0],data[-1][:,1]))]
        #data[-1][:,1] = np.sqrt(data[-1][:,1]+1)
        #print(np.lexsort((data[-1][:,0],data[-1][:,1])))
        #print(data[-1])
    
    
    #fig.set_size_inches(10,5)
    
    # now compute minimum energy phase
    # first need a list of all the fA samples
    
    pt_set = set() # set of all points in phase space to consider
    for j in range(len(data)):
        tmplist = []
        for k in range(data[j].shape[0]):
            tmplist.append(list(data[j][k,0:2]))
        #print(tmplist)
        pt_set = pt_set.union(set(tuple(i) for i in tmplist))
    
    pt_array = np.array(list(pt_set))
    #print(pt_array)
    
    tri = Delaunay(pt_array)
    
    # ax.scatter(pt_array[:,0],pt_array[:,1])
    #plt.show()
    # ax.scatter(mainarray[loc,0], mainarray[loc,1], c=mainarray[loc,2]-2100,s=500.0, cmap="viridis",alpha = 0.5)

    convex_hull = []
    del_list = []
    for j in range(len(data)):
        del_list.append([])
    for j in range(pt_array.shape[0]):
        Emin = 0.
        idxmin = phase_names.index('DIS')
        tmp = np.where(np.logical_and(data[idxmin][:,0]==pt_array[j,0],data[idxmin][:,1]==pt_array[j,1]))[0]
        if (tmp.shape[0] != 0):
            Emin = data[idxmin][tmp[0],2]
        else:
            Emin = np.inf
        for l in range(len(data)):
            tmp = np.where(np.logical_and(data[l][:,0]==pt_array[j,0],data[l][:,1]==pt_array[j,1]))[0]
            if ( tmp.shape[0] != 0 ):
                idx = tmp[0]
                if ( Emin - data[l][idx,2] > tol ):
                    idxmin = l
                    Emin = data[l][idx,2]
                elif ( np.abs( Emin - data[l][idx,2] ) < tol ): # if energy is equal to DIS phase energy, delete it for simplicity
                        del_list[l].append(idx)
        convex_hull.append([pt_array[j,0],pt_array[j,1],idxmin,Emin])
    #print(convex_hull)
    
    del_list[phase_names.index('DIS')] = []
    for j in range(len(del_list)):
        data[j] = np.delete(data[j],del_list[j],axis=0)
    
    
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
            # if counter_keyword==0:
                # ax.scatter(scatter_data[j][:,0],scatter_data[j][:,1],alpha=1.0)

    # if counter_keyword==0:
        # ax.legend(stable_phase_names,bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
        # plt.triplot(pt_array[:,0], pt_array[:,1], tri.simplices,c='k',alpha=0.3)
    

    if draw_boundaries:
        edge_list = []
        phase_boundaries = [];
        phase_boundary_dict = dict()
        boundary_name = []
        for simpidx in range(tri.simplices.shape[0]):
            # generate all pairs of vertices for edges
            tmpedges = list(itertools.combinations(tri.simplices[simpidx,:],2))
            lObjEdges = []
            allbdypts = []
            for edgeidx in range(len(tmpedges)):
                lObjEdges.append(Edge(tmpedges[edgeidx]))
                allbdypts += lObjEdges[-1].getBdyPts()
            #if len(allbdypts)>0:
            #    print(allbdypts)
            if (len(allbdypts)%2 == 0):
                while (len(allbdypts) > 0):
                    len0=len(allbdypts)
                    for l in range(1,len(allbdypts)):
                        #print(bdypts[0],bdypts[l])
                        #print(sorted(bdypts[0][2]),sorted(bdypts[l][2]))
                        if (sorted(allbdypts[0][2]) == sorted(allbdypts[l][2])):
                            phase_boundaries.append([allbdypts[0][0],allbdypts[l][0]])
                            boundary_name.append( [[stable_phase_names[allbdypts[0][2][0]],\
                                              stable_phase_names[allbdypts[0][2][1]]],\
                                             [stable_phase_names[allbdypts[1][2][0]],\
                                              stable_phase_names[allbdypts[1][2][1]]]])
                            del allbdypts[l]
                            del allbdypts[0]
                            break
                    if not (len(allbdypts) < len0):
                        # if we get this far then there were no matches, which means we have a problem
                        print('Error: domain with possible 4 phase intersection point')
                        break
            elif (len(allbdypts) == 3):
            # iterate over all pairs of phases 
            #    #print(allbdypts)
                bdytmp = []
                slopevec =[]
            #    # compute all three boundaries, then find intersection
            #    #print(allbdypts)
                for pidx in range(3): # iterate over all three pairs of phases
                    tmp = []
                    for eidx in range(3): # iterate over all three edges
                        edgetmp = FindEdgeIntersection(tmpedges[eidx],allbdypts[pidx][2])
                        if (len(edgetmp) == 2): # FindEdgeIntersection can return empty
                            tmp.append(edgetmp)
                    #print(tmp)
                    if (len(tmp) == 2):
                        #if the bdy only has one point, then it is just a corner, ignore it, as it will be taken care of in another line
                        bdytmp.append(tmp)
                #print(bdytmp)
                for pidx in range(len(bdytmp)):
                    slopevec.append((bdytmp[pidx][1][1] - bdytmp[pidx][0][1])/(bdytmp[pidx][1][0] - bdytmp[pidx][0][0]))
                if len(bdytmp) == 3:    
                    # strictly all three boundaries should intersect at one point, but they'll probably be a little off. Instead use average of all three individual intersections
                    favg = 0
                    eavg = 0
                    for l in range(len(allbdypts)-1):
                        for m in range(l+1,len(allbdypts)):
                            fcross=-((bdytmp[l][0][1] - slopevec[l]*bdytmp[l][0][0]) - (bdytmp[m][0][1] - slopevec[m]*bdytmp[m][0][0]))/(slopevec[l]-slopevec[m])
                            ecross=bdytmp[l][0][1] + slopevec[l]*(fcross - bdytmp[l][0][0] )
                            favg += fcross
                            eavg += ecross
                    favg /= (len(allbdypts)*(len(allbdypts)-1))/2.
                    eavg /= (len(allbdypts)*(len(allbdypts)-1))/2.
                    #print(favg, eavg)
                    #print(allbdypts)
                    phaseboundary_temp = []
                    for l in range(len(allbdypts)):
                        phase_boundaries.append([allbdypts[l][0],[favg,eavg]])    
                        phaseboundary_temp+=[[stable_phase_names[allbdypts[l][2][0]],stable_phase_names[allbdypts[l][2][1]]]]
                        boundary_name.append(phaseboundary_temp)
                elif len(bdytmp) == 1:
                    phase_boundaries.append(bdytmp[0])
        
        
        
        boundary_name_flat = []
        for j in range(len(boundary_name)):
            for k in range(len(boundary_name[j])):
                boundary_name[j][k].sort()
                boundary_name_flat.append(boundary_name[j][k])
        unique_boundary_set = set(tuple(i) for i in boundary_name_flat)
        unique_boundary = list(list(i) for i in unique_boundary_set)
        countype = np.zeros(len(unique_boundary))
        for j in range(len(phase_boundaries)):
            
            
            
            ax.plot([phase_boundaries[j][0][0],phase_boundaries[j][1][0]],[phase_boundaries[j][0][1],phase_boundaries[j][1][1]],c='k')
            
            # if len(boundary_name[j])>2:
            #     ax.scatter([phase_boundaries[j][0][0],phase_boundaries[j][1][0]],[phase_boundaries[j][0][1],phase_boundaries[j][1][1]],c=color[len(unique_boundary)+1],marker = '^',alpha = 1)
            # else:  
            #     for k in range(len(phase_boundaries[j])):
            #         if countype[unique_boundary.index(boundary_name[j][k])]==0:
            #             ax.scatter(phase_boundaries[j][k][0],phase_boundaries[j][k][1],\
            #                        c=color[unique_boundary.index(boundary_name[j][k])],\
            #                            marker = '^',alpha = 1,label =f'{boundary_name[j][k][0]}-{boundary_name[j][k][1]}' )
                        
            #             countype[unique_boundary.index(boundary_name[j][k])]+=1
                        
            #         else:
            #              ax.scatter(phase_boundaries[j][k][0],phase_boundaries[j][k][1],\
            #             c=color[unique_boundary.index(boundary_name[j][k])],\
            #                 marker = '^',alpha = 1 )

        
    plt.xlim(0.1,0.38)
    # plt.ylim(1.8,2.25)
    
    counter_keyword+=1
# ax.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)

plt.xlabel('$f_A$')
plt.ylabel('$\epsilon$')
plt.savefig(exportname,dpi = 300)

plt.show()
