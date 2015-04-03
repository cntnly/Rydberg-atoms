# -*- coding: utf-8 -*-
"""
Created on Fri Mar 27 12:02:57 2015

@author: Anonymous
"""
from __future__ import division
import numpy as np
from para import *
from constant import *
import pyximport; pyximport.install()
try:
    from radinte_sqrt import *
except ImportError:
    raise Exception('No radinte module found')
def radinte_atom(atomA, atomAp, I =1):
    """
    radinte_atom(atomA, atomAp, I =1)
    Calculate radinte for pair of atoms"""
    return radinte(atomA.E_radinte, atomA.l, atomAp.E_radinte, atomAp.l, I)

try:
    from wigner import *
except ImportError:
    raise Exception('No wigner module found')
    
def A_Stark (lA, jA, mA, lAp, jAp, mAp):
    """
    A_Stark (lA, jA, mA, lAp, jAp, mAp)
    Calculate angular integral for Stark shift, case F aligned with quantization axis
    """
    S = Wigner3j(jA, 1, jAp, -mA, 0, mAp)   
    A = (1-2*((jA + mA + jAp +1/2)%2))*np.sqrt((2*jA + 1)* (2*jAp +1))*Wigner6j(jA, 1, jAp, lAp, 1/2, lA) * np.sqrt((2*lA +1)*(2*lAp +1))* Wigner3j(lA, 1, lAp, 0,0,0)
    return A*S
    
def A_Stark_atom(atomA, atomAp):
    """
    A_Stark_atom(atomA, atomAp)
    Calculate angular integral for Stark shift, case F aligned with quantization axis
    """
    return A_Stark(atomA.l, atomA.j, atomA.m, atomAp.l, atomAp.j, atomAp.m)
    
def R_Int(pairAB, pairABp):
    """
    R_Int(pairAB, pairABp)
    Calculate the radial integral part of the VdW interaction
    """
    return radinte_atom(pairAB.atom1, pairABp.atom1)*radinte_atom(pairAB.atom2, pairABp.atom2)

def A_Integral(lA, jA, mA, lB, jB, mB, lAp, jAp, mAp, lBp, jBp, mBp, theta, epsilon =1e-10):
    """
    A_Int(lA, jA, mA, lB, jB, mB, lAp, jAp, mAp, lBp, jBp, mBp, theta)
    Calculate the angular part of VdW interaction. Needed Wigner3j and Wigner6j
    """
    S = 0.
    if np.abs(2 - 3* np.sin(theta)**2) > epsilon :
        A = Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
        A += Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
        S += - A* (2 - 3* np.sin(theta)**2)*0.5
    if np.abs(np.sin(theta)) > epsilon:
        A = Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
        A += Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
        S += - A* (np.sin(theta)**2)*1.5
    if np.abs(1-3*np.cos(theta)**2) > epsilon:
        A = Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
        S += A*(1- 3* np.cos(theta)**2)
    
    if np.abs(np.sin(2*theta)) > epsilon :
        A = Wigner3j(jA, 1, jAp, -mA, -1, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
        A += -Wigner3j(jA, 1, jAp, -mA, 1, mAp)* Wigner3j(jB, 1, jBp, -mB, 0, mBp)
        A += Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, -1, mBp)
        A += -Wigner3j(jA, 1, jAp, -mA, 0, mAp)* Wigner3j(jB, 1, jBp, -mB, 1, mBp)
        S += -1.5* np.sin(theta)* np.cos(theta)* np.sqrt(2)* A

    A = (1-2*((jA + mA + jAp +1./2)%2))*np.sqrt((2*jA + 1)* (2*jAp +1))*Wigner6j(jA, 1, jAp, lAp, 1/2, lA) * np.sqrt((2*lA +1)*(2*lAp +1))* Wigner3j(lA, 1, lAp, 0,0,0)
    A *= (1-2*((jB + mB + jBp +1./2)%2))*np.sqrt((2*jB + 1)* (2*jBp +1))*Wigner6j(jB, 1, jBp, lBp, 1/2, lB) * np.sqrt((2*lB +1)*(2*lBp +1))* Wigner3j(lB, 1, lBp, 0,0,0)
    S *= A
    return S

def A_Int(pairAB, pairABp):
    """
    A_Int(pairAB, pairABp)
    Call function A_Integral for two pairs
    """
    return A_Integral(pairAB.atom1.l, pairAB.atom1.j, pairAB.atom1.m, pairAB.atom2.l, pairAB.atom2.j, pairAB.atom2.m, pairABp.atom1.l, pairABp.atom1.j, pairABp.atom1.m, pairABp.atom2.l, pairABp.atom2.j, pairABp.atom2.m, theta)


def Search_VdW(N_list, Not_list, test_term, Choice, delta_n, l1, l2, delta_l, m1, m2, delta_m):
    """
    Search_VdW(N_list, Not_list, test_term, Choice, delta_n, l1, l2, delta_l, m1, m2, delta_m)
    """
    up = False
    N_list_temp = []
    for lA in np.arange(max(0, l1-delta_l), l1+delta_l+0.1, 1):
        for lB in np.arange(max(0, l2-delta_l), l2+delta_l +0.1, 1):
            for jA in np.arange(np.abs(lA-0.5), lA+0.6, 1):
                for jB in np.arange(np.abs(lB -0.5), lB +0.6, 1):
                    for mA in np.arange(m1 -delta_m, m1 +delta_m +0.1, 1):
                        if (-jA <= mA) & (mA <= jA):
                            for mB in np.arange(m2-delta_m, m2+delta_m +0.1,1):
                                if (-jB <= mB) & (mB <= jB):
                                    atomA_temp = Ryd_atom(100, lA, jA, mA)
                                    atomB_temp = Ryd_atom(100, lB, jB, mB)
                                    pair_AB_temp = Ryd_pair(atomA_temp, atomB_temp)
                                    for pair in N_list:
                                        A_Int_temp = A_Int(pair_AB_temp,pair)
                                        if A_Int_temp != 0.:
                                            switch_A = 1
                                            iA = 0
                                            while (iA <= delta_n) & (switch_A != 0):
                                                nA = pair.atom1.n + switch_A*iA
                                                if nA > lA:
                                                    atom_A = Ryd_atom(nA, lA, jA, mA)
                                                    switch_B = 1
                                                    iB = 0
                                                    while (iB <= delta_n) & (switch_B != 0):
                                                        nB = pair.atom2.n + switch_B*iB
                                                        if nB > lB:
                                                            atom_B = Ryd_atom(nB, lB, jB, mB)
                                                            pair_AB = Ryd_pair(atom_A, atom_B)
                                                            if  (pair_AB not in Not_list) & (pair_AB not in N_list_temp):
                                                                R_Int_temp = R_Int(pair_AB, pair)
                                                                VColumn_temp = R_Int_temp * A_Int_temp
                                                                if np.abs(VColumn_temp*VColumn_temp*test_term) >= Choice*np.abs(pair_AB.E_Zeeman - pair.E_Zeeman):
                                                                    N_list_temp.append(pair_AB)
                                                                else:
                                                                    switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                                    iB = 0
                                                        if iB == delta_n:
                                                            up = True
                                                            switch_B = -0.5*(switch_B + np.abs(switch_B))
                                                            iB = 0
                                                        iB += 1
                                                        # end search level for B
                                                if iA == delta_n:
                                                    switch_A = -0.5* (switch_A + np.abs(switch_A))
                                                    iA = 0
                                                iA += 1
    if up :
        print('Should increase delta_n')
    return N_list_temp
    

def Search_Stark(pair, Not_list, Choice_F, delta_n_max):
    """
    Search_Stark(pair, Not_list, Choice_F)
    return list of pairs coupled with 'pair' but not in 'Not_list' stronger than '1/Choice_F' 
    """
    N_list_Stark = []
    atom_1, atom_2 = pair.atom1, pair.atom2
    l1, l2 = atom_1.l, atom_2.l
    j1, j2 = atom_1.j, atom_2.j
    m1, m2 = atom_1.m, atom_2.m

    for lA in arange(np.abs(l1 -1), l1+1.1, 2):
        for jA in arange(np.abs(j1-1), j1+1.1, 1):
            if (jA >= np.abs(lA -0.5)) & (jA <= lA +0.5):
                for mA in arange(m1 -1, m1 +1.1, 1):
                    if (-jA <= mA) & (mA <= jA):
                        atomA_temp = Ryd_atom(100, lA, jA, mA)
                        A_Int_temp = A_Stark_atom(atom_1, atomA_temp)
                        if A_Int_temp!=0:
   #                         A_Int_temp *= coef_F
                            #             if A_Int_temp != 0.:
                            switch_A = 1
                            iA = 0
                            while (iA <= delta_n_max) & (switch_A != 0):
                                nA = n1 + switch_A*iA
                                if nA > lA:
                                    atom_A = Ryd_atom(nA, lA, jA, mA)
                                    if ((A_Int_temp*radinte_atom(atom_1, atom_A))**2)*Choice_F> np.abs(atom_1.En - atom_A.En):
                                        pair_temp = Ryd_pair(atom_A, atom_2)                                    
                                        if (pair_temp not in Not_list) & (pair_temp not in N_list_Stark):                                   
                                            N_list_Stark.append(pair_temp)
                                if iA == delta_n_max:
                                    switch_A = -0.5* (switch_A + np.abs(switch_A))
                                    iA = 0
                                iA += 1 
    for lB in arange(np.abs(l1 -1), l1+1.1, 2):
        for jB in arange(np.abs(j1-1), j1+1.1, 1):
            if (jB >= np.abs(lB -0.5)) & (jB <= lB +0.5):
                for mB in arange(m1 -1, m1 +1.1, 1):
                    if (-jB <= mB) & (mB <= jB):
                        atomB_temp = Ryd_atom(100, lB, jB, mB)
                        A_Int_temp = A_Stark_atom(atom_2, atomB_temp)
                        if A_Int_temp!=0:
   #                         A_Int_temp *= coef_F
                            #             if A_Int_temp != 0.:
                            switch_B = 1
                            iB = 0
                            while (iB <= delta_n_max) & (switch_B != 0):
                                nB = n1 + switch_B*iB
                                if nB > lB:
                                    atom_B = Ryd_atom(nB, lB, jB, mB)
         #                           print( (A_Int_temp*radinte_atom(atom_2, atom_B))**2)
         #                           print(np.abs(atom_2.En - atom_B.En))
                                    if ((A_Int_temp*radinte_atom(atom_2, atom_B))**2)*Choice_F> np.abs(atom_2.En - atom_B.En):
                                        pair_temp = Ryd_pair(atom_1,atom_B)
                                        if (pair_temp not in Not_list) & (pair_temp not in N_list_Stark):                                  
                                            N_list_Stark.append(pair_temp)
                                if iB == delta_n_max:
                                    switch_B = -0.5* (switch_B + np.abs(switch_B))
                                    iB = 0
                                iB += 1 
    return N_list_Stark
    
def count_doub(list1, list2):
    """
    Count how many doubled elements in two list
    """
    count =0
    for elm in list1:
        if elm in list2:
            count +=1
    return count