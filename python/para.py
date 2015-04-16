# -*- coding: utf-8 -*-

from __future__ import division
from numpy import pi,cos, sin
import builtins

# Define levels
builtins.n1 = 60
builtins.l1 = n1-1
#s1 = 1/2
#j1 = l1 + 1/2
builtins.m1 = l1

builtins.n2 = 60
builtins.l2 = n2-1     
#s2 =1/2
#j2 = l2 + 1/2
builtins.m2 = l2

builtins.Bfield = 10.e-4 # Magnetic field from experiment 1 Teslta = 10*4 Gauss
builtins.Ffield = 1e-2*sin(0.1*pi/2) # V/cm
builtins.theta = 1*pi/2
builtins.theta_F =1.*pi/2
builtins.phi_F = 0.*pi/2
#theta =  0.95531661812450919

# Setup criterion
delta_n_max = 6
l_max = 2
Choice = 1e7# cut off energy for 1st order term, in Hz
builtins.R_test = 1e-6
Choice2 = 10* Choice

Choice_F = 1e-1 # cut off for Stark shift

if __name__ == '__main__':
    import cal_VdW

