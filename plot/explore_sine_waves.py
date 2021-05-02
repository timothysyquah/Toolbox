#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 14:15:00 2021

@author: tquah
"""





from scipy.fft import fft,ifftn
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

plt.close('all')
# N = 100
# f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, sharex='col', sharey='row')
# xf = np.zeros((N,N))
# xf[0, 3] = 1
# xf[0, N-3] = 1
# Z = ifftn(xf)
# ax1.imshow(xf, cmap=cm.Reds)
# ax4.imshow(np.real(Z), cmap=cm.gray)
# xf = np.zeros((N, N))
# xf[3, 0] = 1
# xf[N-3, 0] = 1
# Z = ifftn(xf)
# ax2.imshow(xf, cmap=cm.Reds)
# ax5.imshow(np.real(Z), cmap=cm.gray)
# xf = np.zeros((N, N))
# xf[5, 10] = 1
# xf[N-5, N-10] = 1
# Z = ifftn(xf)
# ax3.imshow(xf, cmap=cm.Reds)
# ax6.imshow(np.real(Z), cmap=cm.gray)
# plt.show()
plt.figure()
x = np.linspace(0,2*np.pi)
y = np.sin(x)
yy = fft(y)
ytest = ifftn(yy)
plt.plot(x,y,'--')
# plt.plot(x,yy)
# plt.plot(x,ytest)
y = 2*np.sin(x)
yy = fft(y)
ytest = ifftn(yy/2)
plt.plot(x,y)
plt.scatter(x,yy)
plt.scatter(x,ytest)



# import numpy as np
# import matplotlib.pyplot as plt


# def Rotation_2D(x,y,theta_degree):
#     theta = theta_degree*(np.pi/180)
#     xprime = x*np.cos(theta)-y*np.sin(theta)
#     yprime = x*np.sin(theta)+y*np.cos(theta)
#     return xprime, yprime
 # q qq
# x = np.linspace(0,2*np.pi,100)
# y = np.sin(x)
# xprime,yprime = Rotation_2D(x,y,45)
# # plt.plot(x,y)
# plt.plot(xprime,yprime)


# xx,yy = np.meshgrid()
