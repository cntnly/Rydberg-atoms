# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
#import numexpr as ne
import sys
sys.path.append('/home/2kome/Desktop/testdeck/test1/python')

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
coef_F = Ffield*e*e/(4*np.pi*epsilon_0*h)
test_term = coef/(R_test**3)
test_term = test_term*test_term

print(pair_12)
print('theta = {0}'.format(theta*180/pi))
print('B_field = {0} G'.format(Bfield*1e4))


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
                            print atomA_temp
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
                                                        print atomB_temp
                                                        pair_temp = Ryd_pair(atomA_temp, atomB_temp)
                                                        if pair_temp not in N_list:
                                                            N_list.append(pair_temp)
                
print('Degenerated pairs ={0}'.format(len(N_list)))

#Search 1st order coupling terms
N_list1 = Search_VdW(N_list, N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)
#N_list2 = Search_VdW([pair_34],N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)
#N_list3 = Search_VdW([pair_56],N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)
#N_list4 = Search_VdW([pair_78],N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)

def count_doub(list1, list2):
    count =0
    for elm in list1:
        if elm in list2:
            count +=1
    return count
N_list += N_list1
print ('N_list = {0}'.format(len(N_list)))
       

for pair in N_list:
    if pair_invert(pair) not in N_list:
        N_list.append(pair_invert(pair))
print ('1st order terms: {0}'.format(len(N_list)))

    
# Search for 2nd order terms
#test_term2 = coef/(R_test2**3)
#test_term2 = test_term2*test_term2
#N_list2 =  Search_VdW(N_list1, N_list, test_term2, Choice2, delta_n_max2, l1, l2, 2, m1,m2, 2)
#print ('Few Few !!!')
#print ('2st order terms: {0}'.format(len(N_list2)))


#N_list = N_list + N_list2 # sure that elements of N_list are not in N_list2
#for pair in N_list:
#    if pair_invert(pair) not in N_list:
#        print(pair)
#        N_list.append(pair_invert(pair))

# Search for 3nd order terms
#test_term3 = coef/(R_test3**3)
#test_term3 = test_term3*test_term3
#N_list3 =Search_VdW(N_list2, N_list, test_term3, Choice3, delta_n_max2, l1, l2, 3, m1,m2, 3)
#print ('Few Few Few!!!')
#print ('3st order terms: {0}'.format(len(N_list3)))

Union_list = N_list # sure that elements of N_list are not in N_list2
for pair in Union_list:
    if pair_invert(pair) not in Union_list:
        Union_list.append(pair_invert(pair))
        
        
print ('Matrix size: {0}'.format(len(Union_list)))
raw_input("Press enter to continue")
#theta = 0* pi/2
#print(theta*180/pi)


##Create base for testing
#Union_list =[]
#for l1 in [0,1]:
#    for j1 in arange(abs(l1 - 0.5), abs(l1 + 0.6)):
#        for m1 in arange(-j1, j1+1):
#            atom1temp = Ryd_atom(60, l1, j1,m1)
#            for l2 in [0,1]:
#                for j2 in arange(abs(l2 - 0.5), abs(l2 + 0.6)):
#                    for m2 in arange(-j2, j2+1):
#                        atom2temp = Ryd_atom(60, l2,j2,m2)
#                        Union_list.append(Ryd_pair(atom1temp, atom2temp))
                        
# sort list as energy order
Union_list = sorted(Union_list, key = lambda energy: energy.E_Zeeman)
print(len(Union_list))                    
figure(1)
clf()
plot([elm.E_Zeeman - pair_12.E_Zeeman for elm in Union_list])
# Create interaction Matrix
length = len(Union_list)
# create mask
V_total = np.zeros((length,length))
V_R = np.zeros_like(V_total)
V_A = np.zeros_like(V_total)
# mask to calculate only lower triangular matrix
mask1 = np.tril(np.ones_like(V_total))!=0

#vectorize necessary functions
rad_vec = np.vectorize(radinte)
A_vec = np.vectorize(A_Integral)
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

figure(2)
clf()
imshow(V_total,)
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
out_coef = np.empty((4,R_num))
out_vector = np.empty((length, length))
R = np.logspace(log10(R_min), log10(R_max), num = R_num)
for i,elm in enumerate(R):
    #out_egr[i] = np.linalg.eigvalsh(EI + V_total* coef/(elm**3)) - pair_12.E_Zeeman
    out_egr[i] , out_vector = np.linalg.eigh(EI + 1e18*V_total* coef/(elm**3))
    out_egr[i] = out_egr[i] - pair_12.E_Zeeman*1e-9
    out_coef[0,i] = np.argmax(abs(out_vector[index1,:]))
    out_coef[1,i] = np.argmax(abs(out_vector[index2,:]))
#    out_coef[2,i] = np.argmax(abs(out_vector[index3,:]))
#    out_coef[3,i] = np.argmax(abs(out_vector[index4,:]))


out_egr1 = np.asarray([out_egr[i, out_coef[0,i]] for i in range(R_num)])
#out_egr3 = np.asarray([out_egr[i, out_coef[2,i]] for i in range(R_num)])
#out_egr4 = np.asarray([out_egr[i, out_coef[3,i]] for i in range(R_num)])


figure(3)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
semilogx(R, out_egr1)
#semilogx(R, out_egr[:,index3])
#semilogx(R, out_egr[:,index4])

#fit data
def VdW(r,C):
    return C*(r**-6)
from scipy.optimize import curve_fit

#popt,pcov = curve_fit(VdW, R[100:], abs(out_egr1[100:]))
#popt2,pcov = curve_fit(VdW, R[80:], abs(out_egr2[80:]))
#figure(4)
#clf()
#loglog(R, VdW(R, popt), R,abs(out_egr1))
#loglog(R, VdW(R, popt2),R,abs(out_egr2))

#popt,pcov = curve_fit(VdW, R[100:], abs(out_egr[:,index1][100:]))
#popt2,pcov2 = curve_fit(VdW, R[100:], abs(out_egr[:,index2][100:]))
#popt3,pcov3 = curve_fit(VdW, R[100:], abs(out_egr[:,index3][100:]))
#popt4,pcov4 = curve_fit(VdW, R[100:], abs(out_egr[:,index4][100:]))
popt1,pcov1 = curve_fit(VdW, R[100:], out_egr1[100:], p0=100)
figure(4)
clf()
loglog(R, abs(VdW(R, popt1)), R,abs(out_egr1))
print(popt1)

#popt3,pcov3 = curve_fit(VdW, R[100:], abs(out_egr3[100:]- pair_56.E_Zeeman+ pair_12.E_Zeeman))
#popt4,pcov4 = curve_fit(VdW, R[100:], abs(out_egr4[100:]- pair_78.E_Zeeman+ pair_12.E_Zeeman))
if index2 != index1:
    out_egr2 = np.asarray([out_egr[i, out_coef[1,i]] for i in range(R_num)])
    figure(3)
    semilogx(R, out_egr2)
    popt2,pcov2 = curve_fit(VdW, R[100:], out_egr2[100:], p0=100)
    figure(4)
    loglog(R, abs(VdW(R, popt2)), R,abs(out_egr2))
    print(popt2)
    #loglog(R, VdW(R, popt3), R,abs(out_egr[:,index3]))
    #loglog(R, VdW(R, popt4), R,abs(out_egr[:,index4]))


#print(popt2*1e27)
#print(popt3*1e27)
#print(popt4*1e27)

def pow_fit(r,C,n):
    return C*(r**-n)
popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:])
figure(4)
clf()
loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1))
