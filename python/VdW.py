# -*- coding: utf-8 -*-

from __future__ import division
import numpy as np
import numexpr as ne
import os
sys.path.append('D:/NTL/Simulation/Python/python')

from para import *
from constant import * # define constant to be used
# import radinte
import pyximport; pyximport.install()
try:
    from radinte_log import *
except ImportError:
    raise Exception('No radinte module found')
def radinte_pair(atomA, atomAp, I =1):
    """
    radinte_pair(atomA, atomAp, I =1)
    Calculate radinte for pair of atoms"""
    return radinte(atomA.E_radinte, atomA.l, atomAp.E_radinte, atomAp.l, I)

try:
    from wigner import *
except ImportError:
    raise Exception('No wigner module found')
    
def R_Int(pairAB, pairABp):
    """
    R_Int(pairAB, pairABp):
    Calculate the radial integral part of the VdW interaction
    """
    return radinte_pair(pairAB.atom1, pairABp.atom1)*radinte_pair(pairAB.atom2, pairABp.atom2)

def A_Integral(lA, jA, mA, lB, jB, mB, lAp, jAp, mAp, lBp, jBp, mBp, theta):
    """
    A_Int(lA, jA, mA, lB, jB, mB, lAp, jAp, mAp, lBp, jBp, mBp, theta)
    Calculate the angular part of VdW interaction. Needed Wigner3j and Wigner6j
    """
    A = Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
    A += Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
    S = - A* (2 - 3* np.sin(theta)**2)*0.5
    A = Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
    A += Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
    S += - A* (np.sin(theta)**2)*1.5
    A = Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
    S += A*(1- 3* np.cos(theta)**2)

    A = Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
    A = -Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
    A = Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
    A = -Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
    S += 1.5* np.sin(theta)* np.cos(theta)* np.sqrt(2)* A

    A = (1-2*((jA + mA + jAp +1/2)%2))*np.sqrt((2*jA + 1)* (2*jAp +1))*Wigner6j(jA, 1, jAp, lAp, 1/2, lA) * np.sqrt((2*lA +1)*(2*lAp +1))* Wigner3j(lA, 1, lAp, 0,0,0)
    A *= (1-2*((jB + mB + jBp +1/2)%2))*np.sqrt((2*jB + 1)* (2*jBp +1))*Wigner6j(jB, 1, jBp, lBp, 1/2, lB) * np.sqrt((2*lB +1)*(2*lBp +1))* Wigner3j(lB, 1, lBp, 0,0,0)
    S *= A
    return S

def A_Int(pairAB, pairABp):
    """
    A_Int(pairAB, pairABp)
    Call function A_Integral for two pairs
    """
    return A_Integral(pairAB.atom1.l, pairAB.atom1.j, pairAB.atom1.m, pairAB.atom2.l, pairAB.atom2.j, pairAB.atom2.m, pairABp.atom1.l, pairABp.atom1.j, pairABp.atom1.m, pairABp.atom2.l, pairABp.atom2.j, pairABp.atom2.m, theta)


# Calculate Energy
atom_1 = Ryd_atom(n1, l1, j1, m1)
atom_2 = Ryd_atom(n2, l2, j2, m2)
pair_12 = Ryd_pair(atom_1, atom_2)

coef = a_0*a_0*e*e/(4*np.pi*epsilon_0*h)
test_term = coef/(R_test**3)
test_term = test_term*test_term

test_term2 = coef/(R_test2**3)
test_term2 = test_term2*test_term2

N_list = []
N_list2 = []

# Search for 1st order coupling terms
for lA in arange(np.abs(l1 -1), l1+1.1, 2):
    for lB in arange(np.abs(l2 -1), l2+1.1, 2):
        for jA in arange(np.abs(j1-1), j1+1.1, 1):
            if (jA >= np.abs(lA -0.5)) & (jA <= lA +0.5):
                for jB in arange(np.abs(j2-1), j2+1.1, 1):
                    if (jB >= np.abs(lB-0.5)) & (jB<= lB +0.5):
               #         ljA = choose_lj(lA, jA)
               #         ljB = choose_lj(lB, jB)
                        for mA in arange(m1 -1, m1 +1.1, 1):
                            if (-jA <= mA) & (mA <= jA):
                                for mB in arange(m2-1, m2+1.1,1):
                                    if (-jB <= mB) & (mB <= jB):
                                        atomA_temp = Ryd_atom(100, lA, jA, mA)
                                        atomB_temp = Ryd_atom(100, lB, jB, mB)
                                        pair_AB_temp = Ryd_pair(atomA_temp, atomB_temp)
                                        A_Int_temp = A_Int(pair_AB_temp,pair_12)
                                        if A_Int_temp != 0.:
                                            switch_A = 1
                                            iA = 0
                                            while (iA <= delta_n_max) & (switch_A != 0):
                                                nA = n1 + switch_A*iA
                                                if nA > lA:
                                                    atom_A = Ryd_atom(nA, lA, jA, mA)
           #                                         EA = En(nA, ljA)
                                                    # search level for nB
                                                    switch_B = 1
                                                    iB = 0
                                                    while (iB <= delta_n_max) & (switch_B != 0):
                                                        nB = n2 + switch_B*iB
                                                        if nB > lB:
                                                            atom_B = Ryd_atom(nB, lB, jB, mB)
                                                            pair_AB = Ryd_pair(atom_A, atom_B)

                                            # EB = En(nB, ljB)
                                                           # E_AB = EA + EB
                                                            R_Int_temp = R_Int(pair_AB, pair_12)
                                                            VColumn_temp = R_Int_temp * A_Int_temp
                                                            if np.abs(VColumn_temp*VColumn_temp*test_term) >= Choice*np.abs(pair_AB.E_Zeeman - pair_12.E_Zeeman):
                                                                N_list.append(pair_AB)
                                                            else:
                                                                switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                                iB = 0
                                                        if iB == delta_n_max:
                                                            up = True
                                                            switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                            iB = 0
                                                        iB += 1

                                                    # end search level for B
                                                if iA == delta_n_max:
                                                    switch_A = -0.5* (switch_A + np.abs(switch_A))
                                                    iA = 0
                                                iA += 1
print ('Few')
for pair in N_list:
    if pair_invert(pair) not in N_list:
        N_list.append(pair_invert(pair))
print ('1st order terms: {0}'.format(len(N_list)))
if up:
    print('Should increase delta_max')
    
# Search for 2nd order terms
for lA in arange(np.abs(np.abs(l1 -1)-1), l1+2.1, 2):
    for lB in arange(np.abs(np.abs(l2 -1)-1), l2+2.1, 2):
        for jA in arange(np.abs(lA-0.5), lA+0.6, 1):
            for jB in arange(np.abs(lB -0.5), lB +0.6, 1):
                #        ljA = choose_lj(lA, jA)
                #        ljB = choose_lj(lB, jB)
                for mA in arange(m1 -2, m1 +2.1, 1):
                    if (-jA <= mA) & (mA <= jA):
                        for mB in arange(m2-2, m2+2.1,1):
                            if (-jB <= mB) & (mB <= jB):
                                atomA_temp = Ryd_atom(100, lA, jA, mA)
                                atomB_temp = Ryd_atom(100, lB, jB, mB)
                                pair_AB_temp = Ryd_pair(atomA_temp, atomB_temp)
                                for pair in N_list:
                                    A_Int_temp = A_Int(pair_AB_temp,pair)
                                    if A_Int_temp != 0.:
                                        switch_A = 1
                                        iA = 0
                                        while (iA <= delta_n_max) & (switch_A != 0):
                                            nA = pair.atom1.n + switch_A*iA
                                            if nA > lA:
                                                atom_A = Ryd_atom(nA, lA, jA, mA)
                                                #                                         EA = En(nA, ljA)
                                                # search level for nB
                                                switch_B = 1
                                                iB = 0
                                                while (iB <= delta_n_max) & (switch_B != 0):
                                                    nB = pair.atom2.n + switch_B*iB
                                                    if nB > lB:
                                                        atom_B = Ryd_atom(nB, lB, jB, mB)
                                                        pair_AB = Ryd_pair(atom_A, atom_B)
                                                        if  (pair_AB not in N_list) & (pair_AB not in N_list2):
                                                            R_Int_temp = R_Int(pair_AB, pair)
                                                            VColumn_temp = R_Int_temp * A_Int_temp
                                                            if np.abs(VColumn_temp*VColumn_temp*test_term2) >= Choice2*np.abs(pair_AB.E_Zeeman - pair.E_Zeeman):
                                                                N_list2.append(pair_AB)
                                                            else:
                                                                switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                                iB = 0
                                                    if iB == delta_n_max:
                                                        up2 = True
                                                        switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                        iB = 0
                                                    iB += 1
                                                    # end search level for B
                                                if iA == delta_n_max:
                                                    switch_A = -0.5* (switch_A + np.abs(switch_A))
                                                    iA = 0
                                                iA += 1
print ('Few Few !!!')
print ('2st order terms: {0}'.format(len(N_list2)))
if up2:
    print ('Should increase delta_max2')

Union_list = N_list + N_list2 # sure that elements of N_list are not in N_list2
for pair in Union_list:
    if pair_invert(pair) not in Union_list:
        Union_list.append(pair_invert(pair))
print ('Matrix size: {0}'.format(len(Union_list)))

# sort list as energy order
Union_list = sorted(Union_list, key = lambda energy: energy.E_Zeeman)
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


V_A[mask1] = A_vec(Pair1_l1, Pair1_j1, Pair1_m1, Pair1_l2, Pair1_j2, Pair1_m2,Pair2_l1, Pair2_j1, Pair2_m1, Pair2_l2, Pair2_j2, Pair2_m2, theta)
# Choose only elements different from 0
mask2 = V_A[mask1] !=0
V_R[V_A!=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)*rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
V_total[V_A!=0] = V_R[V_A!=0]*V_A[V_A!=0]
V_total = V_total + V_total.T - np.diag(V_total.diagonal())

figure(2)
clf()
imshow(V_total)

# Zero-th Energy
EI = np.diag(np.asarray([elm.E_Zeeman for elm in Union_list]))

R_max = 3.0e-5 #um
R_min = 1.0e-6
#R = R_max
R_num = 100 # step in logarithm
out_egr = np.empty((R_num, length))
R = np.logspace(log10(R_min), log10(R_max), num = R_num)
for i,elm in enumerate(R):
    out_egr[i] = np.linalg.eigvalsh(EI + V_total* coef/(elm**3)) - pair_12.E_Zeeman
index1 = Union_list.index(pair_12)
index2 = Union_list.index(pair_invert(pair_12))
figure(3)
clf()
plot(R, out_egr[:,index1])
plot(R, out_egr[:,index2])


#fit data
def VdW(r,C):
    return C*(r**-6)
from scipy.optimize import curve_fit
popt,pcov = curve_fit(VdW, R[40:], out_egr[:,index1][40:])
loglog(R, VdW(R, popt))

