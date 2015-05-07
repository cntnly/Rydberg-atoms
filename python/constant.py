# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 11:50:49 2015

@author: Anonymous
"""
from __future__ import division
import numpy as np
from scipy.constants import codata
from scipy.constants import h, hbar, c, mu_0, epsilon_0, m_e, e 
k_B = codata.value('Boltzmann constant')
a_0 = codata.value('atomic unit of length') 
mu_B = codata.value('Bohr magneton')
u = codata.value('atomic mass constant')
Ry = codata.value('Rydberg constant')
m_Rb87 = 86.909180520*u
alpha = codata.value('fine-structure constant')
Ry_Rb87 = Ry/(1+m_e/m_Rb87) # Rydberg constant with mass correction

#from para import *

# calculate quantum defect
#delta_0 = np.array([3.1311804, 2.6548849, 2.6416737,1.34809171, 1.34646572, 0.0165192, 0.0165437, 0])
#delta_2 = np.array([0.1784, 0.2900, 0.2950, - 0.60286, - 0.59600, - 0.085, - 0.086, 0])
delta_0, delta_2 = np.zeros(8), np.zeros(8)
def delta(n, lj):
    """
    delta(n.lj)
    Quantum defect for given n and lj for n>20, cf T.Gallager "Rydberg Atoms"
    """
    return delta_0[lj] + delta_2[lj]/((n - delta_0[lj])**2)
def En(n, lj):
    """
    En(n,lj) in Hz
    Energy of Rydberg level with quantum defect
    """
    return - Ry_Rb87*c/(n**2)
def E_radinte(n, lj):
    """
    E_radinte(n, lj) in atomic unit
    energy needed for radinte
    """
    return - 0.5/(1+ m_e/m_Rb87)/(n**2)
def choose_lj(l, j):
    """
    choose_lj(l, j)
    Define function to choose lj    
    """
    if l == 0:
        if j == 1./2:
            return 0
    if l == 1.:
        if j == 1./2:
            return 1
        else:
            return 2
    if l == 2.:
        if j == 3./2:
            return 3
        else:
            return 4
    if l == 3.:
        if j == 5./2:
            return 5
        else:
            return 6
    else:
        return 7


class Ryd_atom(object):
    """
    Ryd_atom(n, l, m)
    Define Rydberg atom with its n, l, m, which ignores fine structure. This holds true for circular series atoms
    Return En
    """
    def __init__(self, n, l, m):
        if l>n:
            raise Exception('l > n')
        if np.abs(m) > l:
            raise Exception('|m| > l')
        self.n = n
        self.l = l
        self.m = m
        self.En = En(n, 7) # 10 to ignore the quantum defect
        self.E_radinte = E_radinte(n, 7)
#        self.E_Zeeman = self.En + Zeemanshift(l, m, Bfield)
#    def E_Zeeman(self, B_field):
#        return self.En + Zeemanshift(self.l, self.m, B_field)
    def __repr__(self):
        if (self.n-1 == self.l)&(self.m ==self.l):
            return 'atom {0:0.0f}C'.format(self.n)
        return "atom {0}, {1}, {2}".format(self.n, self.l, self.m)
    def __eq__(self, other):
        return (self.n, self.l, self.m) == (other.n, other.l, other.m)
        
class Ryd_pair (object):
    def __init__(self, atom1, atom2):
        """
        Define a pair of atom 1 & 2. Return E12
        """
        if isinstance(atom1, Ryd_atom) & isinstance(atom2, Ryd_atom):
            self.atom1 = atom1
            self.atom2 = atom2
        else:
            raise Exception('input should be of Ryd_atom class')
            pass
        self. E = atom1.En + atom2.En
#        self.E_Zeeman = atom1.E_Zeeman + atom2.E_Zeeman
        
            
    def __repr__(self):
        return '{0} and {1}'.format(self.atom1, self.atom2)
    def __eq__(self, other):
        return (self.atom1 == other.atom1) & (self.atom2 == other.atom2)
def pair_invert(pair):
    """ pair_invert(pair)
    Invert pair of atoms"""
    return Ryd_pair(pair.atom2, pair.atom1)
