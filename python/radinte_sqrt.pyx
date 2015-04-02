from __future__ import division
#import numpy as np
#cimport
cimport cython
#cimport numpy as np
from libc.math cimport exp, sqrt, pow
#cython: boundscheck=False
#cython: wraparound=False

#Adapte for Cython
@cython.cdivision(True)
def radinte(double E1, double L1, double E2, double L2, double I, double H = 0.01):
    """
    radinte(double E1, double L1, double E2, double L2, double I, double H = 0.01)
	Calculate Numerov integration for two states <E1,L1| R**I|E2,L2> with square root step H (default 0.01)
	Test formula
	<n,l | r| n, l-1> = 3/2n Sqrt[n^2-l^2]
	"""
	
	
    cdef double AL1 = (2*L1 + 0.5)*(2*L1 + 1.5)
    cdef double AL2 = (2*L2 + 0.5)*(2*L2 + 1.5)
    cdef double EE1 = -2. * E1
    cdef double EE2 = -2. * E2
    cdef double R01 = sqrt(1./EE1)
    cdef double R02 = sqrt(1./EE2)
    R01 = R01*2.*(R01 + 15.)
    R02 = R02*2.*(R02 + 15.)
    
    cdef double C1 = 12./H/H
    cdef double C2 = 2.*C1
    cdef double C3 = 10.
    cdef int I1 = 1   # switch to stop the integral
    cdef int I2 = 1
    cdef double R = R01
    if R01 < R02:
        R = R02
    cdef double U = sqrt(R)    
    cdef double S1 = 0.
    cdef double S2 = 0.
    cdef double S12 = 0.
    cdef double W11 = 1.e-10
    cdef double W21 = 1.e-10
    cdef double W12 = W11*(1. -2.*H*U*EE1)
    cdef double W22 = W21*(1. -2.*H*U*EE2)
    cdef double X1 = pot(U,EE1,AL1)
    cdef double X2 = pot(U,EE2,AL2)
    U = U - H
    cdef double Y1 = pot(R, EE1, AL1)
    cdef double Y2 = pot(R, EE2, AL2)
    
    while (U > 0.01) & ((I1 != 0) | (I2 != 0)):
        U = U - H
        R = U*U
        A = W11*(X1 - C1)
        X1 = Y1
        Y1 = pot(U, EE1, AL1)
        if R <= R01:
            W11 = W12
            W12 = (A + W11*(C3*X1 +C2))/(C1 - Y1)
            if (I1 ==1) & (X1 <= 0.) & (Y1 > 0.):
                I1 =2
            if I1 == 2:
                if W11/W12 < 1.:
                    I1 = 0
                    W11 = 0.
                    W12 = 0.
        A = W21 * (X2 -C1)
        X2 = Y2
        Y2 = pot(U, EE2, AL2)
        if R <= R02:
            W21 = W22
            W22 = (A + W21*(C3*X2 +C2))/(C1 - Y2)
            if (I2 == 1) & (X2 <= 0.) & (Y2 > 0.):
                I2 = 2
            if I2 == 2 :
                if W21/W22 < 1.:
                    I2 = 0
                    W21 = 0.
                    W22 = 0.
        S1 += (W12*U)*(W12*U)
        S2 += (W22*U)*(W22*U)
        S12+= W12*W22*(U**(2.*I+2.))
    return S12/sqrt(S1)/sqrt(S2)
cdef double pot(double R, double E, double A):
    """Define the potential as pure Coulombic. Polarization and Fine structure terms can be added"""
    R = R*R
    return -8. + 4.*E*R + A/R



