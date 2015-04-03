# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
#import numexpr as ne
import sys
sys.path.append('C:/Users/r14/Documents/GitHub/test/python')

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
    


# Calculate Energy
atom_1 = Ryd_atom(n1, l1, j1, m1)
atom_2 = Ryd_atom(n2, l2, j2, m2)
pair_12 = Ryd_pair(atom_1, atom_2)
pair_34 =Ryd_pair(Ryd_atom(n1,l1,j1,-m1),Ryd_atom(n2, l2, j2, -m2))
pair_56 =Ryd_pair(Ryd_atom(n1,l1,j1,m1),Ryd_atom(n2, l2, j2, -m2))
pair_78 =Ryd_pair(Ryd_atom(n1,l1,j1,-m1),Ryd_atom(n2, l2, j2, m2))


coef = a_0*a_0*e*e/(4*np.pi*epsilon_0*h)
test_term = coef/(R_test**3)
test_term = test_term*test_term

print(pair_12)
print('theta = {0} deg'.format(theta*180/pi))
print('B_field = {0} G'.format(Bfield*1e4))
print('F_field = {0} V/cm'.format(Ffield))
Ffield = Ffield*100 # V/m
coef_F = Ffield*e*a_0/h

#check degeneracy
N_list = [pair_12]
for lA in np.arange(max(l1-2, 0), l1+2.1,1):
    if lA < n1:
        for jA in np.arange(max(j1- 2, 0.5), j1 + 2.1,1):
            if (jA >= np.abs( lA - 0.5)) & (jA <= lA + 0.5):
                for mA in np.arange(m1-2, m1+ 2.1,1):
                    if np.abs(mA) <= jA:                 
                        try:
                            atomA_temp = Ryd_atom(n1, lA, jA, mA)
                        except Exception:
                            print(lA, jA, mA)
                        if np.abs(atomA_temp.E_Zeeman - atom_1.E_Zeeman) < 1e-10:
                            #print (atomA_temp)
                            for lB in np.arange(max(l2-2, 0), l2+2.1):
                                if lB < n2:
                                    for jB in np.arange(max(j2- 2, 0.5), j2 + 2.1):
                                        if (jB >= np.abs(lB - 0.5)) & (jB <= lB + 0.5):
                                            for mB in np.arange(m2-2, m2+ 2.1):
                                                if np.abs(mB) <= jB:
                                                    try:
                                                        atomB_temp = Ryd_atom(n2, lB, jB, mB)
                                                    except Exception:
                                                        print (lB, jB, mB)
                                                    if np.abs(atomB_temp.E_Zeeman - atom_2.E_Zeeman) < 1e-10:
                                                        #print (atomB_temp)
                                                        pair_temp = Ryd_pair(atomA_temp, atomB_temp)
                                                        if pair_temp not in N_list:
                                                            N_list.append(pair_temp)
deg=len(N_list)
print('Degenerated pairs ={0}'.format(len(N_list)))

#Search 1st order coupling terms
N_list1 = Search_VdW(N_list, N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)


N_list += N_list1
print ('N_list = {0}'.format(len(N_list)))
       

for pair in N_list:
    if pair_invert(pair) not in N_list:
        N_list.append(pair_invert(pair))
print ('1st order terms: {0}'.format(len(N_list)))


# Search for temp coupled by Stark effect
N_list_Stark=[]
Choice_F = (coef_F*coef_F)/Choice_F
for elm in N_list:
    N_list_Stark_temp = Search_Stark(elm, N_list + N_list_Stark, Choice_F, delta_n_max)
    N_list_Stark += N_list_Stark_temp

print('1st order Stark terms: {0}'.format(len(N_list_Stark)))

count_doub(N_list, N_list_Stark)

#============= Create base 
#for n in arange(atom_1.n -15, atom_1.n +16, 1):
#    for l in arange(max(0, atom_1.l -3), min(atom_1.l +3.1,n),1):
#        for j in arange(abs(l-0.5), l+0.6,1):
#            for m in arange(-j, j+0.1,1):
#                if abs(m -atom_1.m) <3:
#                    try:
#                        atom_temp = Ryd_atom(n,l,j,m)
#                    except:
#                        print(n, l,j,m)                            
#                    if abs(atom_temp.En - atom_1.En) < 140e9:
#                        if  atom_temp not in N_list:
#                            N_list.append(atom_temp)
#                            
#                            
#Union_list = N_list
#=============================

Union_list = N_list + N_list_Stark # sure that elements of N_list are not in N_list2
for pair in Union_list:
    if pair_invert(pair) not in Union_list:
        Union_list.append(pair_invert(pair))
        
        
print ('Matrix size: {0}'.format(len(Union_list)))
#try:
#    raw_input("Press enter to continue")
#except:
#    input("Press Enter to continue")
                        
# sort list as energy order
Union_list = sorted(Union_list, key = lambda energy: energy.E_Zeeman)
print(len(Union_list))                    
figure(1)
clf()
subplot(2,2,1)
plot([(elm.E_Zeeman - pair_12.E_Zeeman)*1e-9 for elm in Union_list],'-o')
xlabel('List N°')
ylabel('Rel. energy (GHz)')
# Create interaction Matrix
length = len(Union_list)
# create mask
V_total = np.zeros((length,length))
V_R = np.zeros_like(V_total)
V_A = np.zeros_like(V_total)
V_Stark1 = np.zeros_like(V_total)
V_Stark2 = np.zeros_like(V_total)
# mask to calculate only lower triangular matrix
mask1 = np.tril(np.ones_like(V_total))!=0

#vectorize necessary functions
rad_vec = np.vectorize(radinte)
A_vec = np.vectorize(A_Integral)
A_Stark_vec = np.vectorize(A_Stark)

V_row = np.asarray([(elm.atom1.n, elm.atom1.l, elm.atom1.j, elm.atom1.m, elm.atom1.E_radinte, elm.atom2.n, elm.atom2.l, elm.atom2.j, elm.atom2.m, elm.atom2.E_radinte) for elm in Union_list])
Pair1_n1, Pair2_n1 = meshgrid(V_row[:,0],V_row[:,0])
Pair1_l1, Pair2_l1 = meshgrid(V_row[:,1],V_row[:,1])
Pair1_j1, Pair2_j1 = meshgrid(V_row[:,2],V_row[:,2])
Pair1_m1, Pair2_m1 = meshgrid(V_row[:,3],V_row[:,3])
Pair1_Erad1, Pair2_Erad1 = meshgrid(V_row[:,4],V_row[:,4])
Pair1_n2, Pair2_n2 = meshgrid(V_row[:,5],V_row[:,5])
Pair1_l2, Pair2_l2 = meshgrid(V_row[:,6],V_row[:,6])
Pair1_j2, Pair2_j2 = meshgrid(V_row[:,7],V_row[:,7])
Pair1_m2, Pair2_m2 = meshgrid(V_row[:,8],V_row[:,8])
Pair1_Erad2, Pair2_Erad2 = meshgrid(V_row[:,9],V_row[:,9])

# Saving memory
Pair1_n1, Pair1_l1, Pair1_j1, Pair1_m1, Pair1_Erad1 = Pair1_n1[mask1], Pair1_l1[mask1], Pair1_j1[mask1], Pair1_m1[mask1], Pair1_Erad1[mask1]
Pair1_n2, Pair1_l2, Pair1_j2, Pair1_m2, Pair1_Erad2 = Pair1_n2[mask1], Pair1_l2[mask1], Pair1_j2[mask1], Pair1_m2[mask1], Pair1_Erad2[mask1]
Pair2_n1, Pair2_l1, Pair2_j1, Pair2_m1, Pair2_Erad1 = Pair2_n1[mask1], Pair2_l1[mask1], Pair2_j1[mask1], Pair2_m1[mask1], Pair2_Erad1[mask1]
Pair2_n2, Pair2_l2, Pair2_j2, Pair2_m2, Pair2_Erad2 = Pair2_n2[mask1], Pair2_l2[mask1], Pair2_j2[mask1], Pair2_m2[mask1], Pair2_Erad2[mask1]

V_A[mask1] = A_vec(Pair1_l1, Pair1_j1, Pair1_m1, Pair1_l2, Pair1_j2, Pair1_m2, Pair2_l1, Pair2_j1, Pair2_m1, Pair2_l2, Pair2_j2, Pair2_m2, theta)
# Choose only elements different from 0
mask2 = V_A[mask1] !=0
V_R[V_A!=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)*rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
V_total[V_A!=0] = V_R[V_A!=0]*V_A[V_A!=0]
V_total = V_total + V_total.T - np.diag(V_total.diagonal())
V_total = V_total*1e-9 # convert to GHz.m^3

figure(1)
subplot(2,2,2)
imshow(V_total,)
colorbar()

V_R = np.zeros_like(V_total)
V_A = np.zeros_like(V_total)
V_A[mask1] = A_Stark_vec(Pair1_l1, Pair1_j1, Pair1_m1, Pair2_l1, Pair2_j1, Pair2_m1)*(Pair1_n2 == Pair2_n2) *(Pair1_l2 == Pair2_l2)*(Pair1_j2 == Pair2_j2)*(Pair1_m2 == Pair2_m2)
mask2 = V_A[mask1] !=0
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)
    V_Stark1[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark1 = V_Stark1 + V_Stark1.T - np.diag(V_Stark1.diagonal())
    V_Stark1 = V_Stark1*1e-9
subplot(2,2,3)
imshow(V_Stark1)
colorbar()

V_R = np.zeros_like(V_total)
V_A = np.zeros_like(V_total)
V_A[mask1] = A_Stark_vec(Pair1_l2, Pair1_j2, Pair1_m2, Pair2_l2, Pair2_j2, Pair2_m2)*(Pair1_n1 == Pair2_n1) *(Pair1_l1 == Pair2_l1)*(Pair1_j1 == Pair2_j1)*(Pair1_m1 == Pair2_m1)
mask2 = V_A[mask1] !=0
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
    V_Stark2[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark2 = V_Stark2 + V_Stark2.T - np.diag(V_Stark2.diagonal())
    V_Stark2 = V_Stark2*1e-9
subplot(2,2,4)
imshow(V_Stark2)
colorbar()

# Zero-th Energy
EI = np.diag(np.asarray([elm.E_Zeeman for elm in Union_list]))*1e-9

R_max = 30 #in um
R_min = 1.
#R = R_max
R_num = 200 # step in logarithm
index1 = Union_list.index(pair_12)
if pair_invert(pair_12) != pair_12:
    index2 = Union_list.index(pair_invert(pair_12))
else:
    index2 = index1
#index2 = Union_list.index(pair_invert(pair_34))
#index3 = Union_list.index(pair_invert(pair_56))
#index4 = Union_list.index(pair_invert(pair_78))

out_egr = np.empty((R_num, length))
out_coef = np.empty((2,R_num))
out_vector = np.empty((R_num,length, length))
R = np.logspace(log10(R_min), log10(R_max), num = R_num)
for i,elm in enumerate(R):
    #out_egr[i] = np.linalg.eigvalsh(EI + V_total* coef/(elm**3)) - pair_12.E_Zeeman
    out_egr[i] , out_vector[i] = np.linalg.eigh(EI + 1e18*V_total* coef/(elm**3) + coef_F*(V_Stark1 + V_Stark2))
    out_egr[i] = out_egr[i] - pair_12.E_Zeeman*1e-9
    out_coef[0,i] = np.argmax(abs(out_vector[i][index1,:]))
    out_coef[1,i] = np.argmax(abs(out_vector[i][index2,:]))

#    idx = index1
#    coef_idx = np.where(np.argmax(abs(out_vector[i]), axis=0)== idx)
#    if len(coef_idx[0]) !=0:
#        out_coef[0,i] = coef_idx[0][np.argmin(abs(coef_idx[0]-idx))]
#    if len(coef_idx[0]) == 0:
#        while len(coef_idx[0]) == 0:
#            idx -= 1
#            coef_idx = np.where(np.argmax(abs(out_vector[i]), axis=0)== idx)
#        out_coef[0,i] = coef_idx[0][0]
#    idx = index2
#    coef_idx = np.where(np.argmax(abs(out_vector[i]), axis=0)== idx)
#    if len(coef_idx[0]) !=0:
#        out_coef[1,i] = coef_idx[0][np.argmin(abs(coef_idx[0]-idx))]
#    if len(coef_idx[0]) == 0:
#        while len(coef_idx[0]) == 0:
#            idx -= 1
#            coef_idx = np.where(np.argmax(abs(out_vector[i]), axis=0)== idx)
#        out_coef[1,i] = coef_idx[0][0]
        
#if deg >1:
#    out_egr1 = out_egr[:,out_coef[0,R_num-1]]
#    out_egr2 = out_egr[:,out_coef[1,R_num-1]]
#else:
#    out_egr1 = out_egr[:, index1]
#    out_egr2 = out_egr[:, index2]
    
out_egr1 = [out_egr[i, out_coef[0,i]] for i in range(R_num)]
out_egr2 = [out_egr[i, out_coef[1,i]] for i in range(R_num)]



figure(2)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
loglog(R, asarray(out_egr), R, asarray(out_egr1), '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')
#semilogx(R, out_egr[:,index3])
#semilogx(R, out_egr[:,index4])
figure(3)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
semilogx(R, out_egr, R, out_egr1, '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')

#fit data
def VdW(r,C):
    return C*(r**-6)
def pow_fit(r,C,n):
    return C*(r**-n)

#popt,pcov = curve_fit(VdW, R[100:], abs(out_egr1[100:]))
#popt2,pcov = curve_fit(VdW, R[80:], abs(out_egr2[80:]))
#figure(4)
#clf()
#loglog(R, VdW(R, popt), R,abs(out_egr1))
#loglog(R, VdW(R, popt2),R,abs(out_egr2))

from scipy.optimize import curve_fit

offset, off_vec = np.linalg.eigh(EI + coef_F*(V_Stark1 + V_Stark2))
offset = offset - pair_12.E_Zeeman*1e-9
offset1 = offset[np.argmax(abs(off_vec[index1,:]))]
offset2 = offset[np.argmax(abs(off_vec[index2,:]))]

popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:]-offset1, p0=(100,6))
figure(4)
clf()
loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1-offset1), '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')
print('power fit = {0}'.format(popt1))
print('R_1 MHz = {0} um'.format(R[np.max(np.where(abs(out_egr1-offset1)>1e-3))]))


#popt3,pcov3 = curve_fit(VdW, R[100:], abs(out_egr3[100:]- pair_56.E_Zeeman+ pair_12.E_Zeeman))
#popt4,pcov4 = curve_fit(VdW, R[100:], abs(out_egr4[100:]- pair_78.E_Zeeman+ pair_12.E_Zeeman))
if index2 != index1:
  #  out_egr2 = np.asarray([out_egr[i, out_coef[1,i]] for i in range(R_num)])
    #out_egr2 = out_egr[:, index2]
    figure(3)
    semilogx(R, out_egr2, 'wo')
    popt2,pcov2 = curve_fit(pow_fit, R[100:], out_egr2[100:], p0=(100,6))
    figure(4)
    loglog(R, abs(pow_fit(R, *popt2)), R,abs(out_egr2-offset2),'wo')
    print(popt2)
  





#popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:])
#figure(5)
#clf()
#loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1))
#print(popt1)