#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:01:17 2021

@author: tquah
"""
import numpy as np
import matplotlib.pyplot as plt

def multiple_formatter(denominator=2, number=np.pi, latex='\pi'):
    def gcd(a, b):
        while b:
            a, b = b, a%b
        return a
    def _multiple_formatter(x, pos):
        den = denominator
        num = np.int(np.rint(den*x/number))
        com = gcd(num,den)
        (num,den) = (int(num/com),int(den/com))
        if den==1:
            if num==0:
                return r'$0$'
            # if num==1:
            #     return r'$%s$'%latex
            # elif num==-1:
            #     return r'$-%s$'%latex
            else:
                return r'$%s\times %s$'%(num,latex)
        else:
            if num==1:
                return r'$\frac{%s}{%s}$'%(latex,den)
            elif num==-1:
                return r'$\frac{-%s}{%s}$'%(latex,den)
            else:
                return r'$\frac{%s%s}{%s}$'%(num,latex,den)
    return _multiple_formatter


class Zbounds():
    def __init__(self,C,Nsc):
        self.C = C
        self.Nsc = Nsc
    def z1(self):
        return self.C*np.power(self.Nsc,-1/2)
    def z22(self):
        return self.C
    def z21(self):
        return self.C*self.C
    def z31(self):
        return np.sqrt(self.C)
Nsc = np.linspace(5,50,100)

C = np.linspace(0,3,100)
constant = 1#6*np.sqrt(6)
    
plt.close('all')
plt.figure()
ax = plt.gca()
C1 = np.linspace(0,1,1000)
C2 = np.linspace(1,2,1000)
Nsc_array = np.arange(5,45,5,dtype=int)
Cplot1 = C1*constant
Cplot2 = C2*constant



zbound1 = Zbounds(C1,0)
z21_ = zbound1.z21()
z31_ = zbound1.z31()
plt.plot(Cplot1,z21_,'k')
plt.plot(Cplot1,z31_,'k' )



zbound2 = Zbounds(C2,0)
z22_ = zbound2.z22()
plt.plot(Cplot2,z22_,'k' )

z = np.linspace(-10,10,100)

# plt.plot(np.ones_like(z),z,'--k' )





alpha0 = 1.0
for i in range(0,len(Nsc_array)):
    Cfull = np.linspace(1/np.sqrt(Nsc_array[i]),2,2000)
    Cfullplot = Cfull*constant
    zfull = Zbounds(Cfull,Nsc_array[i])
    
    z1_ = zfull.z1()
    plt.plot(Cfullplot,z1_,'k',alpha = alpha0-0.75/len(Nsc_array)*i )
plt.xlabel(r'$b^3/v_0$')
# ax.xaxis.set_major_locator(plt.MultipleLocator(6*np.sqrt(6)))
# # ax.yaxis.set_minor_locator(plt.MultipleLocator(3*np.sqrt(6)))
# ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter(number =6*np.sqrt(6), latex = ' 6 \sqrt{6}')))

plt.ylabel('$z$')


plt.ylim(-0.1,2.1)
plt.text(1.0*constant,1.5,'Fully Stretched DB',fontsize=15,ha='center', va='center',color='r')
plt.text(1.5*constant,1.0,'LB',fontsize=15,ha='center', va='center',color='r')
plt.text(0.3*constant,0.3,'DB',fontsize=15,ha='center', va='center',color='r')
plt.text(1.5*constant,0.0,'DC',fontsize=15,ha='center', va='center',color='r')
plt.arrow(1.75*constant,0.85,0,-0.7,head_width=0.05*constant,head_length=0.15,length_includes_head=True,color = 'b')
plt.text(1.85*constant,0.68,'$N_{sc}$',fontsize=15,ha='center', va='center',color='b')

plt.tight_layout()

plt.savefig('/home/tquah/Figures/zvsC.png',dpi = 300)


# zstar = 1/np.sqrt(Nsc)
# z1 = 1/Nsc
# plt.close('all')

# plt.figure()
# plt.plot(zstar,Nsc,'--k',label = '$z_1$ - Loosely grafted combs')
# plt.plot(z1,Nsc,'-.k',label = '$z^{\star}$ - Densely grafted combs')

# plt.xlabel('$z$')
# plt.ylabel('$N_{sc}$')
# plt.xlim(0.05,0.5)
# plt.ylim(5,50)
# plt.legend()

# plt.figure()
# cn = np.linspace(0,2,100)

# Nsc_list = np.arange(5,50,5)
# alpha_list = 0.5+(0.5/len(Nsc_list))


# for i in range(len(Nsc_list)):
#     zstar = cn*Nsc_list[i]**(-1/2)

#     alphaval = 0.5+(0.5/len(Nsc_list))*i

#     if i==len(Nsc_list)-1:    
#         plt.plot(cn,zstar,'-r',label = 'z*',alpha= alphaval)
#     else:
#         plt.plot(cn,zstar,'-r',alpha= alphaval)

        
# gzero = np.where(cn<=1)[0]





# zstarstar = np.ones_like(cn)*cn
# zstarstar[gzero] = cn[gzero]*cn[gzero]
# #zstarstar2 = cn
# zstarstarstar=np.sqrt(cn)




# # plt.plot(cn,zstar,'.r',label = '$z^{\star}$ - Loosely grafted combs')
# plt.plot(cn,zstarstar,'-.g',label = '$z^{\star \star}$ - Densely grafted combs')
# plt.plot(cn,zstarstarstar,'-.b',label = '$z^{\star \star \star}$ - Brush')



# plt.xlabel('$c/\sqrt{N}$')
# plt.ylabel('$z$')
# plt.legend()
# plt.savefig('/home/tquah/Figures/zvsc.pdf',dpi = 300)

# plt.figure()
# cn = np.linspace(0,2,100)

# Nsc_list = np.arange(5,50,5)
# alpha_list = 0.5+(0.5/len(Nsc_list))


# for i in range(len(Nsc_list)):
#     z0 = cn*(Nsc_list[i]**(-1/2))
#     z1 = cn*(Nsc_list[i]**(-1.0))

#     alphaval = 0.5+(0.5/len(Nsc_list))*i

#     if i==len(Nsc_list)-1:    
#         plt.plot(cn,z0,'-r',label = 'LC-DC',alpha= alphaval)
#         plt.plot(cn,z1,'-k',label = 'DC-LB',alpha= alphaval)

#     else:
#         plt.plot(cn,z0,'-r',alpha= alphaval)
#         plt.plot(cn,z1,'-k',alpha= alphaval)

        





# zstarstar = np.ones_like(cn)*cn
# zstarstar[gzero] = cn[gzero]*cn[gzero]
# #zstarstar2 = cn
# zstarstarstar=cn




# # plt.plot(cn,zstar,'.r',label = '$z^{\star}$ - Loosely grafted combs')
# # plt.plot(cn,zstarstar,'-.g',label = '$z^{\star \star}$ - Densely grafted combs')
# plt.plot(cn,zstarstarstar,'-.b',label = 'X-Brush')



# plt.xlabel('$c/\sqrt{N}$')
# plt.ylabel('$z$')
# plt.legend()
# plt.savefig('/home/tquah/Figures/zvsc.pdf',dpi = 300)

# # plt.figure()


# # Nbb = 40
# # Nsc = 20
# # b = 1
# # l= 1
# # v =1 
 
# # z = np.linspace(0,2,100)

# # c = (l*z*Nsc)**(3/2)/v**(1/2)/(1+z*Nsc)/(1+z*Nsc)*Nbb
# # plt.plot(z,c)

# # c = z**(5/2)*Nsc**(3/2)/(1+z*Nsc)/(1+z*Nsc)*Nbb
# # plt.plot(z,c)

# # c = v*z**(3)*Nsc**(3/2)/(1+z*Nsc)/(1+z*Nsc)/(b*l)**3/2*Nbb
# # plt.plot(z,c)


