# -*- coding: utf-8 -*-

from __future__ import division
from numpy import pi
# Define levels
n1 = 60
l1 = 59
s1 = 1/2
j1 = l1 + 1/2
m1 = j1
n2 = 60
l2 = 59     
s2 =1/2
j2 = l2 + 1/2
m2 = j2

Bfield = 0.001*8.17e-4 # Magnetic field from experiment 1 Teslta = 10*4 Gauss
Ffield = 0*0.1 # V/cm
theta = 0*pi/2
#theta =  0.95531661812450919

# Setup criterion
delta_n_max = 6
l_max = 2
Choice = 1e06# cut off energy for 1st order term, in Hz
R_test = 1e-6
Choice2 = 10* Choice
Choice3 = 10000* Choice2
R_test2 = 1* R_test
R_test3 = R_test2
delta_n_max2 = 2

Choice_F = 1e0 # cut off for Stark shift
