# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 00:01:29 2015

@author: 2kome
"""

from __future__ import division
import numpy as np
#cimport
#cimport cython
#cimport numpy as np
#from libc.math cimport exp, sqrt, pow
#cython: boundscheck=False
#cython: wraparound=False

#Adapte for Cython
#@cython.cdivision(True)
def wave_func( n1,  L1, delta=0,  H = 0.01):
    """
    radinte( E1,  L1,  E2,  L2,  I,  H = 0.01)
	calculate Numerov integration for two states <E1,L1| R**I|E2,L2> with square root step H (default 0.01)
	radinte(E1, L1, E2, L2, I, H = 0.01)    """
    radius = []
    wave = []
    E1 = -0.5/(n1-delta)**2
 #   E2 = -0.5/n2/n2
    AL1 = (2*L1 + 0.5)*(2*L1 + 1.5)
 #   AL2 = (2*L2 + 0.5)*(2*L2 + 1.5)
    EE1 = -2. * E1
#    EE2 = -2. * E2
    R01 = np.sqrt(1./EE1)
#    R02 = np.sqrt(1./EE2)
    R01 = R01*2.*(R01 + 15.)
 #   R02 = R02*2.*(R02 + 15.)
    
    C1 = 12./H/H
    C2 = 2.*C1
    C3 = 10.
    I1 = 1   # switch to stop the integral
#    I2 = 1
    R = R01
#    if R01 < R02:
 #     R = R02
    U = np.sqrt(R)    
    S1 = 0.
 #   S2 = 0.
 #   S12 = 0.
    W11 = 1.e-10
#    W21 = 1.e-10
    W12 = W11*(1. -2.*H*U*EE1)
#    W22 = W21*(1. -2.*H*U*EE2)
    X1 = pot(U,EE1,AL1)
 #   X2 = pot(U,EE2,AL2)
    U = U - H
    Y1 = pot(R, EE1, AL1)
 #   Y2 = pot(R, EE2, AL2)
    
    while (U > H) & (I1 != 0):#| (I2 != 0)):
        U = U - H
        R = U*U
        A = W11*(X1 - C1)
        X1 = Y1
        Y1 = pot(U, EE1, AL1)
        if R <= R01:
            W11 = W12
            W12 = (A + W11*(C3*X1 +C2))/(C1 - Y1)
            radius.append(R)
            wave.append(W12/(R**0.75))
            if (I1 ==1) & (X1 <= 0.) & (Y1 > 0.):
                I1 =2
            if I1 == 2:
                if W11/W12 < 1.:
                    I1 = 0
                    W11 = 0.
                    W12 = 0.
   #     A = W21 * (X2 -C1)
   #     X2 = Y2
   #     Y2 = pot(U, EE2, AL2)
   #     if R <= R02:
   #         W21 = W22
   #         W22 = (A + W21*(C3*X2 +C2))/(C1 - Y2)
   #         if (I2 == 1) & (X2 <= 0.) & (Y2 > 0.):
   #             I2 = 2
   #         if I2 == 2 :
   #             if W21/W22 < 1.:
   #                 I2 = 0
   #                 W21 = 0.
   #                 W22 = 0.
        S1 += (W12*U)*(W12*U)
   #     S2 += (W22*U)*(W22*U)
   #     S12+= W12*W22*pow(U,(2.*I+2.))
    return asarray(radius), asarray(wave)/np.sqrt(S1*2*H)
def  pot( R,  E,  A):
    """Define the potential as pure Coulombic. Polarization and Fine structure terms can be added"""
    R = R*R
    return -8. + 4.*E*R + A/R



