# -*- coding: utf-8 -*-

from __future__ import division
# Define levels
n1 = 50
l1 = 49
s1 = 1/2
j1 = l1 + 1/2
m1 = j1
n2 = 51
l2 = 49     
s2 =1/2
j2 = l2 + 1/2
m2 = j2

Bfield = 8.17e-4 # Magnetic field from experiment

theta = 0


# Setup criterion
delta_n_max = 6
l_max = 2
Choice = 1e5 # cut off energy for 1st order term, in Hz
R_test = 1e-6
Choice2 = 1* Choice
R_test2 = 2* R_test
delta_n_max2 = 4
