# -*- coding: utf-8 -*-

from __future__ import division
from numpy import pi
# Define levels
n1 = 60
l1 = n1-1
#s1 = 1/2
#j1 = l1 + 1/2
m1 = l1

n2 = 60
l2 = n2-1     
#s2 =1/2
#j2 = l2 + 1/2
m2 = l2

Bfield = 10.e-4 # Magnetic field from experiment 1 Teslta = 10*4 Gauss
Ffield = 1e-1 # V/cm
theta = 1*pi/2
theta_F = 0.001*pi/2
phi_F = 0.001*pi/2
#theta =  0.95531661812450919

# Setup criterion
delta_n_max = 6
l_max = 2
Choice = 1e7# cut off energy for 1st order term, in Hz
R_test = 1e-6
Choice2 = 10* Choice

Choice_F = 1e-1 # cut off for Stark shift
