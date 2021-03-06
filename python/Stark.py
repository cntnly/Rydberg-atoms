# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 16:31:05 2015

@author: Anonymous
"""

from __future__ import division
from numpy import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt

#import numexpr as ne
import sys
sys.path.append('C:/Users/r14/Documents/GitHub/test/python')
sys.path.append('D:/NTL/Simulation/Python/python')

from imp import reload

import para
from para import *

reload(para)
from para import *

import constant
from constant import * # define constant to be used

reload(constant)
from constant import *
# import radinte
import pyximport; pyximport.install()
try:
    from radinte_sqrt import *
except ImportError:
    raise Exception('No radinte module found')
try:
    from functions import *
except ImportError:
    raise Exception('No functions module found')
coef_F = e*a_0/h
Choice_F = (coef_F*coef_F)/Choice_F
atom = Ryd_atom(n1, l1, m1)
print(atom)
print('B_field = {0} G'.format(Bfield*1e4))
#print('theta_F = {0} deg'.format(theta_F*180/pi))
#print('phi_F = {0} deg'.format(phi_F*180/pi))
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(0)
ax = fig.gca(projection='3d')
ax.plot([0,0],[0,0],[0,1],'red')
ax.plot([0,np.sin(theta)],[0,0],[0,np.cos(theta)],'-o')
ax.plot([0,np.sin(theta_F)*np.cos(phi_F)],[0,np.sin(theta_F)*np.sin(phi_F)],[0,np.cos(theta_F)])
plt.show()

N_list = [atom]
def Search_Stark_level(atom, Not_list, Choice, delta_n_max):
    n, l, m = atom.n,atom.l, atom.m
    N_list = []
    for lA in arange(np.abs(l -1), l+1.1, 1):
#        for jA in arange(np.abs(j-1), j+1.1, 1):
#            if (jA >= np.abs(lA -0.5)) & (jA <= lA +0.5):
                for mA in arange(m -1, m +1.1, 1):
                    if (-lA <= mA) & (mA <= lA):
                        atomA_temp= Ryd_atom(100, lA, mA)
                        A_Int_temp = A_Stark_atom(atom, atomA_temp)
                        if A_Int_temp!=0:
   #                         A_Int_temp *= coef_F
                            #             if A_Int_temp != 0.:
                            switch_A = 1
                            iA = 0
                            while (iA <= delta_n_max) & (switch_A != 0):
                                nA = n + switch_A*iA
                                if nA > lA:
                                    atom_A = Ryd_atom(nA, lA, mA)
                                    if atom_A not in (N_list + Not_list):
                                        if ((A_Int_temp*radinte_atom(atom, atom_A))**2)*Choice_F> np.abs(atom.En - atom_A.En):
                                            N_list.append(atom_A)
                                if iA == delta_n_max:
                                   switch_A = -0.5* (switch_A + np.abs(switch_A))
                                   iA = 0
                                iA += 1 
    return N_list
#N_list1= Search_Stark_level(atom, [], Choice_F, delta_n_max)
#print('N_list1 = {0}'.format(len(N_list1)))
#count_doub(N_list, N_list1)
#N_list += N_list1
#N_list2 = []
#for Atom in N_list:
#    N_list_temp = Search_Stark_level(Atom, N_list + N_list2, Choice_F, delta_n_max)
#    N_list2 += N_list_temp
#Union_list = N_list + N_list2

#============= Create base 
def stark_base(N_list, delta_n, delta_l, delta_m, delta_E):
    for n in arange(atom.n -delta_n, atom.n +delta_n+0.1, 1):
        for l in arange(max(0, atom.l -delta_l), min(n,atom.l +delta_l +.1),1):
            for m in arange(-l, l+0.1,1):
                if abs(m -atom.m) < delta_m:
                    atom_temp = Ryd_atom(n, l, m)
                    if abs(atom_temp.En - atom.En) < delta_E:
                        if  atom_temp not in N_list:
                            N_list.append(atom_temp)
    return N_list

Union_list = stark_base(N_list, 3,30,5,1e11)
#=============================

Union_list = sorted(Union_list, key = lambda energy: energy.En)
print('Union_list= {0}'.format(len(Union_list)) )
figure(1)
clf()
plot([(elm.En - atom.En)*1e-9 for elm in Union_list],'-o')
xlabel('List N°')
ylabel('Rel. energy (GHz)')
# Create interaction Matrix
length = len(Union_list)
# create mask
V_Stark = np.zeros((length,length))
V_R = np.zeros_like(V_Stark)
V_A = np.zeros_like(V_Stark)

rad_vec = np.vectorize(radinte)
A_Stark_vec = np.vectorize(A_Stark)

V_row = np.asarray([(elm.n, elm.l, elm.m, elm.E_radinte) for elm in Union_list])
Atom1_n, Atom2_n = meshgrid(V_row[:,0],V_row[:,0])
Atom1_l, Atom2_l = meshgrid(V_row[:,1],V_row[:,1])
#Atom1_j, Atom2_j = meshgrid(V_row[:,2],V_row[:,2])
Atom1_m, Atom2_m = meshgrid(V_row[:,2],V_row[:,2])
Atom1_Erad, Atom2_Erad = meshgrid(V_row[:,3],V_row[:,3])
# mask to calculate only lower triangular matrix
mask1 = np.tril(np.ones_like(V_Stark))!=0
Atom1_n, Atom1_l, Atom1_m, Atom1_Erad = Atom1_n[mask1], Atom1_l[mask1], Atom1_m[mask1], Atom1_Erad[mask1]
Atom2_n, Atom2_l, Atom2_m, Atom2_Erad = Atom2_n[mask1], Atom2_l[mask1], Atom2_m[mask1], Atom2_Erad[mask1]
V_A[mask1] = A_Stark_vec(Atom1_l, Atom1_m, Atom2_l, Atom2_m)
mask2 = V_A[mask1] !=0
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Atom1_Erad[mask2], Atom1_l[mask2], Atom2_Erad[mask2], Atom2_l[mask2], 1)
    V_Stark[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark = V_Stark + V_Stark.T - np.diag(V_Stark.diagonal())
    V_Stark = V_Stark*1e-6
figure(2)
clf()
imshow(np.real(V_Stark))
colorbar()

F_max = 2 #in V/cm
F_min = 0.
F_num = 50 # step in linear 
index = Union_list.index(atom)

# Zero-th Energy
EI = np.diag(np.asarray([elm.E_Zeeman for elm in Union_list]))*1e-6


out_egr = np.empty((F_num, length))
out_coef = np.empty(F_num)
out_coef2 = np.empty(F_num)
out_vector = np.empty((F_num,length, length))
F = np.linspace(F_min, F_max, num = F_num)
for i,elm in enumerate(F*100):
    #out_egr[i] = np.linalg.eigvalsh(EI + V_total* coef/(elm**3)) - pair_12.E_Zeeman
    out_egr[i] , out_vector[i] = np.linalg.eigh(EI + elm*coef_F*V_Stark )
    out_egr[i] = out_egr[i] - atom.En*1e-6
    out_coef[i] = np.argmax((out_vector[i][index,:])**2)
 #   out_coef2[i] = np.argmax((out_vector[i][19,:])**2)
out_egr1 = np.asarray([out_egr[i, out_coef[i]] for i in range(F_num)])
#out_egr2 = np.asarray([out_egr[i, out_coef2[i]] for i in range(F_num)])

#out_egr1 = out_egr[:, index]

figure(3)
clf()
plot(F, out_egr1, '+', F, out_egr,)# F, out_egr2,'wo')
#ylim(-500.,1000.)
xlabel('F (V/cm)')
ylabel('Rel. energy (MHz)')


def fit_fun(F, A,n, offset):
    return A*F**2 + offset
from scipy.optimize import curve_fit
popt1,pcov1 = curve_fit(fit_fun, F[:100], out_egr1[:100], p0=(1, 2, 0))

figure(4)
clf()
plot(F, asarray([out_vector[i, index, out_coef[i]] for i in range(F_num)])**2)
ylim(0,1.1)

figure(5)
clf()
plot(F, out_egr1, 'wo', F, fit_fun(F, *popt1))
xlabel('F (V/cm))')
ylabel('Rel. energy (MHz)')
#title('Stark shift for {atm}, Coef = {coef}'.format(atm =atom, coef = popt1[0]))
text(0.8, 0, '{atm} \n $\Delta E = {A:10.2f}\\times F^{{\;{n:10.2f}}} $ MHz/(V/cm)${{}}^{{{n:0.2f}}}$'.format(atm =atom, A = popt1[0], n = popt1[1]),  fontsize=14, verticalalignment='center')
print('Coef = {a} MHz/(V/cm)**{b}'.format(a= popt1[0], b = popt1[1]))

#def sandwich(vec1, vec2, ope, epsilon =1e-6):
#    c1, c2 = meshgrid(vec1, vec2)
#    c= c1*c2
#    u1, u2 = where(abs(c)> epsilon)
#    S=0
#    for i in range(len(u1)):
#     #   print(c[u1[i],u2[i]],Union_list[u1[i]], Union_list[u2[i]])
#        S += c[u1[i],u2[i]]*ope(Union_list[u1[i]], Union_list[u2[i]])  
#    return S
#dipole = zeros(F_num)
#for i in range(F_num):
# #   u = np.argmax((out_vector[i][19,:])**2)
#    dipole[i]=sandwich(out_vector[i,:,out_coef[i]],out_vector[0,:,out_coef2[i]],lambda atom1, atom2: radinte_atom(atom1,atom2)*A_Stark_atom(atom1,atom2, 0,0), 1e-3)
##dipole2 = zeros(F_num)
##for i in range(F_num):
## #   u = np.argmax((out_vector[i][19,:])**2)
##    dipole2[i]=sandwich(out_vector[i,:,out_coef[i]],out_vector[0,:,out_coef2[i]],lambda atom1, atom2: radinte_atom(atom1,atom2)*A_Stark_atom(atom1,atom2, pi/2,0))
#    
#figure(6);
#plot(F, 1e-6*dipole**2,'-',label='{0}'.format(Union_list[19]))
##plot(F, 1e-6*dipole2**2,'b--')
#legend(loc='best',)
##ylim(-0.1,3.5)
