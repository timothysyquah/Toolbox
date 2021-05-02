#!/usr/bin/env python3

import numpy as np
import sys
import os
sys.path.append("/home/tquah/toolbox_github/lib/")
import glob
#mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
#libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
#sys.path.append(libpath)
import matplotlib.pyplot as plt
import iotools as io
from skimage import filters,measure
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks,peak_prominences
from scipy.stats import cauchy
from skspatial.objects import Plane,Line
from copy import deepcopy
from scipy.optimize import curve_fit
import viztools
import pandas as pd

'''
Method taken/stolen from 
https://stackoverflow.com/questions/52502896/how-can-i-fit-a-good-lorentzian-on-python-using-scipy-optimize-curve-fit


'''


def lorentzian( x, x0, a, gam ):
    return a * gam**2 / ( gam**2 + ( x - x0 )**2)

def multi_lorentz( x, params ):
    off = params[0]
    paramsRest = params[1:]
    assert not ( len( paramsRest ) % 3 )
    return off + sum( [ lorentzian( x, *paramsRest[ i : i+3 ] ) for i in range( 0, len( paramsRest ), 3 ) ] )

def res_multi_lorentz( params, xData, yData ):
    diff = [ multi_lorentz( x, params ) - y for x, y in zip( xData, yData ) ]
    return diff




def load_fields_3D(file,field_id):
    coords, fields = io.ReadDatFile(infile)
    nfields = fields.shape[-1]
    if field_id>nfields:
        print('Warning No Field!')
    else:
        singlefield = fields[:,:,:,field_id]
    return coords,fields,singlefield

def filter_single_field(single_field,resolution):
    edge_sobel = filters.sobel(single_field)
    model = gaussian_kde(np.ravel(edge_sobel))
    x = np.linspace(0,1.0,resolution)
    y = model.pdf(x)
    peaks, _ = find_peaks(y)
    prom,left_bases, right_bases = peak_prominences(y,peaks)
    cutoff = x[right_bases[-2]]
    loc = np.where(edge_sobel>cutoff)
    filtered_image = np.zeros_like(edge_sobel)
    filtered_image[loc] = 1
    labeled_image = measure.label(filtered_image)
    return filtered_image,labeled_image
    

def Plane_fitting(Limage,tol=100):
    unique_planes = np.unique(Limage)
    numpoints = []
    planes = []
    for i in range(1,len(unique_planes)):
        loc = np.where(Limage==unique_planes[i])
        if len(loc[0])>tol:
            numpoints.append(len(loc[0]))
            extractedpoints = np.vstack(loc).transpose()
            plane = Plane.best_fit(extractedpoints)
            planes.append(plane)
    return planes,numpoints
    
def Planes_Check(plane_list,npoints,tol= 0.1):   
    compile_vector_list = [np.array(i.normal) for i in plane_list]      
    vector_array = np.vstack(compile_vector_list)
    if len(plane_list)>2:
        parallel_check = [np.divide(np.array(i.normal),vector_array) for i in plane_list]
        parallel_array = np.vstack(parallel_check)
        meanvalue = np.mean(parallel_array,axis = 1)
        logic = parallel_array-np.vstack(([meanvalue]*3)).transpose()
        diff = np.where(np.abs(logic)>tol/0.5)  
        if len(diff[0]>1):
            print('Error Planes are not Parallel')
            return 0
        else:
            return 1
    else:
        print('Error not enough planes')
        return 0
    

def Reorient_Planes(plane_list,npoints):
    #using plane with most orient the rest
    i= np.where(npoints == np.max(npoints))[0][0]
    signs = np.sign(plane_list[i].normal)
    orient_line = Line(plane_list[i].normal,plane_list[i].point)
    newplane = []
    for i in range(len(plane_list)):
        normalvector = np.abs(np.array(plane_list[i].normal))*signs
        points = np.array(plane_list[i].intersect_line(orient_line))
        newplane.append(Plane(points,normalvector))
    return newplane




def ReorderPlanes(plane_list):
    normalvector = plane_list[0].normal
    point0 = plane_list[0].point
    pointref = 1e6*np.array(normalvector)+point0
    points = np.vstack([i.point for i in plane_list])
    point0 = np.vstack([pointref for i in plane_list])
    
    distance = np.sqrt(np.sum(np.power(np.vstack(points-point0),2),axis = 1))
    index = np.argsort(distance)
    return index
    
def get_distance(plane_list,index):
    dappend = []
    for j in range(1,len(plane_list)-1):
        i = index[j]
        im1 = index[j-1]
        normalline = Line(plane_list[im1].point,direction = plane_list[0].normal) 
        # print(plane_list[i-1].point)
        p0 = np.array(plane_list[im1].point)
        intersection_point = np.array(plane_list[i].intersect_line(normalline))
        distance = np.sqrt(np.sum((intersection_point-p0)**2))
        dappend.append(distance)
    return dappend





def get_angle(x1,x2):
    lhs = np.dot(x1,x2)/np.abs(np.linalg.norm(x1))/(np.linalg.norm(x2))
    return np.arccos(lhs)

def get_average_angle(listofplanes):
    xvector = np.array([0,1,0])
    angle = []
    for i in range(1,len(listofplanes)-1):
        raw_angle = get_angle(listofplanes[i].normal,xvector)*180/np.pi
        angle.append(raw_angle)
    return angle
        


def Replicate_8(coord,singlefield):
    #replacates 8 boxes 2x2x2
    dim = np.shape(coord)[-1]
    
    maxvalue = []
    for i in range(dim):
        maxvalue.append(np.max(coord[:,:,:,i]))
        
        
    singlefield_augment = np.block([[[singlefield,singlefield],\
                                     [singlefield,singlefield]],\
                                    [[singlefield,singlefield],\
                                     [singlefield,singlefield]]])

    coord000 = deepcopy(coord)
    
    coord010 = deepcopy(coord)
    coord010[:,:,:,1] = coord010[:,:,:,1]+maxvalue[1]

    coord100 = deepcopy(coord)
    coord100[:,:,:,0] = coord010[:,:,:,0]+maxvalue[0]
    
    coord110 = deepcopy(coord)
    coord110[:,:,:,1] = coord110[:,:,:,1]+maxvalue[1]
    coord110[:,:,:,0] = coord110[:,:,:,0]+maxvalue[0]

    coord001 = deepcopy(coord)
    coord001[:,:,:,2] = coord001[:,:,:,2]+maxvalue[2]



    coord011 = deepcopy(coord)
    coord011[:,:,:,1] = coord011[:,:,:,1]+maxvalue[1]
    coord011[:,:,:,2] = coord011[:,:,:,2]+maxvalue[2]

    coord101 = deepcopy(coord)
    coord101[:,:,:,0] = coord101[:,:,:,0]+maxvalue[0]
    coord101[:,:,:,2] = coord101[:,:,:,2]+maxvalue[2]

    coord111 = deepcopy(coord)
    coord111[:,:,:,1] = coord111[:,:,:,1]+maxvalue[1]
    coord111[:,:,:,0] = coord111[:,:,:,0]+maxvalue[0]
    coord111[:,:,:,2] = coord111[:,:,:,2]+maxvalue[2]



    coord_augment1 = np.block([[[[coord000]],[[coord010]]],\
                               [[[coord100]],[[coord110]]]])
    coord_augment2 = np.block([[[[coord001]],[[coord011]]],\
               [[[coord101]],[[coord111]]]])
        
    coord_augment = np.concatenate((coord_augment1,coord_augment2),axis = 2)
    return coord_augment,singlefield_augment
    
def Cauchy_Distribution(x,a,b):
    return cauchy.pdf(x,a,b)


def SKAA_Domain(array_in):
    
    xspace = np.linspace(np.min(array_in[:,0]),np.max(array_in[:,0]),1000)
    warning = 0
    try:
        var,pcov = curve_fit(Cauchy_Distribution,array_in[:,0],array_in[:,1],p0=[1.0,0.1])
        domain_space = 2*np.pi/var[0]
        yspace = Cauchy_Distribution(xspace,var[0],var[1])
    except Exception:
        xspace = np.nan
        yspace = np.nan
        domain_space= np.nan
        warning = 1
        pass
    return xspace,yspace,domain_space,warning

def SKAA_Domain_Extra(array_in):
    
    xspace = np.linspace(np.min(array_in[:,0]),np.max(array_in[:,0]),1000)
    warning = 0
    var= [0,0]
    try:
        var,pcov = curve_fit(Cauchy_Distribution,array_in[:,0],array_in[:,1],p0=[1.0,0.1])
        domain_space = 2*np.pi/var[0]
        yspace = Cauchy_Distribution(xspace,var[0],var[1])
    except Exception:
        xspace = np.nan
        yspace = np.nan
        domain_space= np.nan
        warning = 1
        pass
    return xspace,yspace,domain_space,warning,var[0]

def Create_CStress_Tensor(df):
    L = df['L']
    H = df['Hamiltonian.Real']
    list_of_stress = ['StressXX.Real','StressYY.Real','StressZZ.Real',\
                  'StressXY.Real','StressYZ.Real','StressXZ.Real']
    stress_tensor_loc = [[[0,0]],[[1,1]],[[2,2]],\
                     [[0,1],[1,0]],[[1,2],[2,1]],[[0,2],[2,0]]]    
    cauchy_tensor = np.zeros((3,3))
    for i in range(len(list_of_stress)):
        for j in range(len(stress_tensor_loc[i])):
            cauchy_tensor[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress[i]]
    eigenval,eigenvect = np.linalg.eig(cauchy_tensor)

    return L,H,cauchy_tensor,eigenval,eigenvect

def Analyze_Principle_ST(eigenvalues,relative_tolerance):
    loc = [[0,1],[0,2],[1,2]]
    check = np.zeros(3)
    for i in range(len(loc)):
        check[i] = np.isclose(eigenvalues[loc[i]][0],eigenvalues[loc[i]][1],rtol = relative_tolerance)
    return check

def Check_Stress_Operator(infile):
    data_dictionary = dict()
    df = pd.read_csv(infile,delimiter=" ")
    shape = df.shape
    # primary_
    # print('L Eigenvalue_xx Eigenvalue_yy Eigenvalue_zz STATUS')
    datalist = []
    for i in range(shape[0]):
        L, H, cauchy_tensor, eigenval, eigenvect = Create_CStress_Tensor(df.iloc[i])
        check = Analyze_Principle_ST(eigenval,0.05)
        value = len(np.nonzero(check)[0])
        if value==1:
            Status = 'axi-symmetric'
        elif value==3:
            Status = 'isotropic'
        else:
            Status = 'None'
        datalist.append([L,eigenval[0],eigenval[1],eigenval[2],Status])
    return datalist

plt.close('all')
path = '/home/tquah/Figures/'
a = 2
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/0.01_Amplitude_REFINE')
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/SCFT_into_CL/CL/CL_0.01/REFINE/')

list_of_path = glob.glob('L_*')



#
# list_of_path.pop(4)

# list_of_path.pop(5)


# list_of_path = ['L_98','L_100','L_102','L_104']

# list_of_path = ['L_9.2','L_9.22','L_9.23','L_9.24']

outfile = open('Domain_Stats.dat','w+')
outfile.write("L D0_sf D0_pm D0_std nvect_x nvect_y nvect_z angle angle_std \n")
# Check_Stress_Operator('operators_result.dat')


cmap = plt.get_cmap('gnuplot')
colors = [cmap(i) for i in np.linspace(0, 1, 7)]
Domain_Information = []
fig = plt.figure()
ax = plt.gca()


for i in range(len(list_of_path)):
    print(list_of_path[i])
    L = float(list_of_path[i][2:])
    infile = list_of_path[i]+'/DensityOperator.dat'
    coords,fields,singlefield = load_fields_3D(infile,0 )
    aug_coords, aug_single_field = Replicate_8(coords,singlefield)
    filtered_image,labeled_image = filter_single_field(aug_single_field,1000)
    del aug_single_field
    ExportArray = np.stack([filtered_image,labeled_image],axis= 3)
    
    del filtered_image
    viztools.writeVTK(f'filter_L_{L}.vtk', aug_coords, ExportArray)
    
    del ExportArray,aug_coords
    planes,numpoints = Plane_fitting(labeled_image)
    del labeled_image
    newplane = Reorient_Planes(planes,np.array(numpoints))
    
    # passfail = Planes_Check(newplane,numpoints)
    # if passfail==0:
    #     pass
    index = ReorderPlanes(newplane)
    # for i in index:
    #     print(newplane[i].point)
    distance_list = get_distance(newplane,index)
    print(distance_list)
    mean = np.mean(distance_list)
    std = np.std(distance_list)
    period_mu = mean*a*2
    period_std = std*a*2
    array = np.loadtxt(list_of_path[i]+'/SKAA.dat')
    loc = np.where(array[:,0]<2)[0]
    xspace,yspace,domain_space,warning = SKAA_Domain(array[loc,:])
    # if warning==0:
    #     ax.set_yscale('log')
    #     ax.plot(xspace,yspace,color= colors[i])
    #     ax.scatter(array[loc,0],array[loc,1],color = colors[i])
    angle =  get_average_angle(newplane)
    angle_average = np.mean(angle)
    angle_std = np.std(angle)
    Domain_Information.append([L,domain_space,period_mu,period_std,angle_average,angle_std])
    
    outfile.write(f'{L} {domain_space} {period_mu} {period_std} {newplane[0].normal[0]} {newplane[0].normal[1]} {newplane[0].normal[2]} {angle_average} {angle_std} \n')
outfile.close()  
ax.set_ylim(-3,2)  
plt.savefig('SKAA.png',dpi = 300)











# cmap = plt.get_cmap('gnuplot')
# colors = [cmap(i) for i in np.linspace(0, 0.8, 7)]
# Domain_Information = []
# fig = plt.figure()
# ax = plt.gca()
# for i in range(len(list_of_path)):
#     print(list_of_path[i])
#     array = np.loadtxt(list_of_path[i]+'/SKAA.dat')
#     loc = np.where(array[:,0]<2)[0]
#     ax.scatter(array[loc,0],array[loc,1],color = colors[i],label='L='+list_of_path[i][2:])
#     ax.plot(array[loc,0],array[loc,1],color = colors[i])

# ax.set_yscale('symlog')
# ax.set_xlabel('$q(l)$')
# ax.set_xlabel('$S_{AA}(q)/CN$')

# # ax.set_xlim(0,0.4)
# plt.legend()

# plt.savefig('/home/tquah/Figures/SKAA1.png',dpi = 300)
# fig = plt.figure()
# ax = plt.gca()
# for i in range(len(list_of_path)):
#     print(list_of_path[i])
#     array = np.loadtxt(list_of_path[i]+'/SKAA.dat')
#     loc = np.where(array[:,0]<2)[0]
#     xspace,yspace,domain_space,warning,mvalue = SKAA_Domain_Extra(array[loc,:])
#     print(mvalue)
#     if warning==0 and mvalue>0:
#         ax.scatter(array[loc,0],array[loc,1],color = colors[i],label='L='+list_of_path[i][2:])
#         ax.plot(np.ones_like(yspace)*mvalue,yspace,color = colors[i],linestyle = '--',)
#         ax.plot(array[loc,0],array[loc,1],color = colors[i])

# ax.set_yscale('symlog')
# plt.legend()
# ax.set_xlabel('$q(l)$')
# ax.set_xlabel('$S_{AA}(q)/CN$')

# ax.set_xlim(0,0.4)

# plt.savefig('/home/tquah/Figures/SKAA2.png',dpi = 300)
# fig = plt.figure()
# ax = plt.gca()

# for i in range(len(list_of_path)):
#     print(list_of_path[i])
#     array = np.loadtxt(list_of_path[i]+'/SKAA.dat')
#     loc = np.where(array[:,0]<2)[0]
#     xspace,yspace,domain_space,warning = SKAA_Domain(array[loc,:])
#     if warning==0:
#         ax.scatter(array[loc,0],array[loc,1],color = colors[i],label='L='+list_of_path[i][2:])
#         ax.plot(array[loc,0],array[loc,1],color = colors[i])

#         ax.plot(xspace,yspace,color= colors[i],linestyle = '--',alpha = 0.5)
# ax.set_xlim(0,1)

# ax.set_yscale('symlog')
# plt.legend()
# ax.set_xlabel('$q(l)$')
# ax.set_xlabel('$S_{AA}(q)/CN$')

# plt.savefig('/home/tquah/Figures/SKAA3.png',dpi = 300)