#!/usr/bin/env python3

import numpy as np
import pdb
import sys
import os
sys.path.append("/home/tquah/toolbox_github/lib/")
from scipy.fftpack import fftn, fftshift,ifftn
import glob
#mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
#libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
#sys.path.append(libpath)
import matplotlib.pyplot as plt
import iotools as io
from domaintools import DomainAnalyzer
from skimage import feature,measure
from scipy import ndimage as ndi
from skimage import filters
from skimage.data import camera
from skimage.util import compare_images
from scipy.stats import gaussian_kde
from scipy.signal import peak_widths
from scipy.interpolate import CubicSpline
from scipy.signal import chirp, find_peaks, peak_widths,peak_prominences
from scipy.stats import linregress,cauchy
from mpl_toolkits.mplot3d import axes3d
from scipy.linalg import lstsq
from skspatial.objects import Points, Plane,Line
from skspatial.plotting import plot_3d
from copy import deepcopy
from scipy.optimize import curve_fit

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
    

def Plane_fitting(Limage):
    unique_planes = np.unique(Limage)
    numpoints = []
    planes = []
    for i in range(1,len(unique_planes)):
        loc = np.where(Limage==unique_planes[i])
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
    for j in range(1,len(plane_list)):
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
    lhs = np.dot(x1,x2)/np.abs(np.norm(x1))/(np.norm(x2))
    return np.arccos(lhs)



'''
xData, yData = np.loadtxt('HEMAT_1.dat', unpack=True )
yData = yData / max(yData)

generalWidth = 1

yDataLoc = yData
startValues = [ max( yData ) ]
counter = 0

while max( yDataLoc ) - min( yDataLoc ) > .1:
    counter += 1
    if counter > 20: ### max 20 peak...emergency break to avoid infinite loop
        break
    minP = np.argmin( yDataLoc )
    minY = yData[ minP ]
    x0 = xData[ minP ]
    startValues += [ x0, minY - max( yDataLoc ), generalWidth ]
    popt, ier = leastsq( res_multi_lorentz, startValues, args=( xData, yData ) )
    yDataLoc = [ y - multi_lorentz( x, popt ) for x,y in zip( xData, yData ) ]

print popt
testData = [ multi_lorentz(x, popt ) for x in xData ]

'''


def Test_Lorentzian(xData,yData,datax):
    yData = yData / max(yData)
    
    generalWidth = 1e-5
    yDataLoc = yData
    startValues = [ max( yData ) ]
    counter = 0
    
    while max( yDataLoc ) - min( yDataLoc ) > .1:
        counter += 1
        if counter > 20: ### max 20 peak...emergency break to avoid infinite loop
            break
        minP = np.argmin( yDataLoc )
        minY = yData[ minP ]
        x0 = xData[ minP ]
        startValues += [ x0, minY - max( yDataLoc ), generalWidth ]
        popt, ier = leastsq( res_multi_lorentz, startValues, args=( xData, yData ) )
        yDataLoc = [ y - multi_lorentz( x, popt ) for x,y in zip( xData, yData ) ]
    
    return[ multi_lorentz(x, popt ) for x in xData ]
    

def Cauchy_Distribution(x,a,b):
    return cauchy.pdf(x,a,b)


def SKAA_Domain(array_in):
    xspace = np.linspace(np.min(array_in[:,0]),np.max(array_in[:,0]),1000)
    var,pcov = curve_fit(Cauchy_Distribution,array_in[:,0],array_in[:,1],p0=[1.0,0.1])
    domain_space = 2*np.pi/var[0]
    yspace = Cauchy_Distribution(xspace,var[0],var[1])
    return xspace,yspace,domain_space



plt.close('all')



path = '/home/tquah/Figures/'



infile = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/extrarefine/c_9.2_1000/SKAA.dat'





# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/30-1-Averaged_DISLAM/REFINE')
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/0.01_Amplitude_REFINE')
a = 2
list_of_path = glob.glob('L_*')

list_of_path = ['L_6.4','L_8.2','L_9.2','L_9.4','L_10.0']

# for i in range(len(list_of_path)):
#     print(list_of_path[i])

#     infile = list_of_path[i]+'/SKAA.dat'
#     array = np.loadtxt(infile)
#     loc = np.where(array[:,0]<2)[0]
#     xspace = np.linspace(np.min(array[loc,0]),np.max(array[loc,0]),1000)
#     validation = Test_Lorentzian(array[loc,0],array[loc,1],xspace)    
    
#     plt.figure()

#     plt.scatter(array[loc,0],array[loc,1]/np.max(array[loc,1]),label = 'L '+list_of_path[i][2:])
#     plt.plot(array[loc,0],validation,label = 'L '+list_of_path[i][2:])

#     plt.legend()
#     break



Domain_Information = []
plt.figure()
for i in range(len(list_of_path)):
    print(list_of_path[i])
    L = float(list_of_path[i][2:])
    infile = list_of_path[i]+'/DensityOperator.dat'
    array = np.loadtxt(list_of_path[i]+'/SKAA.dat')
    loc = np.where(array[:,0]<2)[0]
    xspace,yspace,domain_space = SKAA_Domain(array[loc,:])
    plt.plot(xspace,yspace)
    coords,fields,singlefield = load_fields_3D(infile,0 )
    filtered_image,labeled_image = filter_single_field(singlefield,1000)
    planes,numpoints = Plane_fitting(labeled_image)
    newplane = Reorient_Planes(planes,np.array(numpoints))
    
    passfail = Planes_Check(newplane,numpoints)
    if passfail==0:
        pass

    index = ReorderPlanes(newplane)
    # for i in index:
    #     print(newplane[i].point)
    
    distance_list = get_distance(newplane,index)
    print(distance_list)
    mean = np.mean(distance_list)
    std = np.std(distance_list)
    
    period_mu = mean*a*2
    period_std = std*a*2
    Domain_Information.append([L,period_mu,period_std,domain_space])

    

