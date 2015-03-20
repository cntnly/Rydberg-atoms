# -*- coding: utf-8 -*-

from __future__ import division
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

Bfield = 0*8.17e-4 # Magnetic field from experiment

theta = 3.1415/2


# Setup criterion
delta_n_max = 16
l_max = 2
Choice = 1e5 # cut off energy for 1st order term, in Hz
R_test = 1e-6
Choice2 = 1* Choice
R_test2 = 2* R_test
delta_n_max2 = 14
