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

from para import *

# calculate quantum defect
delta_0 = np.array([3.1311804, 2.6548849, 2.6416737,
1.34809171, 1.34646572, 0.0165192, 0.0165437, 0])
delta_2 = np.array([0.1784, 0.2900, 0.2950, - 0.60286, - 0.59600, - 0.085, - 0.086, 0])
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
    return - Ry_Rb87*c/((n-delta(n, lj))**2)
def E_radinte(n, lj):
    """
    E_radinte(n, lj) in atomic unit
    energy needed for radinte
    """
    return - 0.5/(1+ m_e/m_Rb87)/((n-delta(n, lj))**2)
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

def g(l, j):
    """
    g(l, j)
    Define Lande g-factor 
    """
    return 1.5 + (3/4 - l*(l+1))/(2*j*(j +1))
def Zeemanshift(l, j, m, Bfield):
    """
    Zeemanshift(l, j, m, Bfield))
    """
    return Bfield * g(l, j)* m * mu_B/h
class Ryd_atom(object):
    """
    Define Rydberg atom with its n, l, j ,m
    Return En, Zeeman shift, E_radinte,"""
    def __init__(self, n, l, j, m):
        self.n = n
        self.l = l
        self. j = j
        self.m = m
        self.lj = choose_lj(l,j)
        self.En = En(n, self.lj)
        self.E_radinte = E_radinte(n, self.lj)
        self.Zeemanshift = Zeemanshift(l, j, m, Bfield)
    def __repr__(self):
        return "atom {0}, {1}, {2} ,{3}".format(self.n, self.l, self.j, self.m)
    def __eq__(self, other):
        return (self.n, self.l, self.j, self.m) == (other.n, other.l, other.j, other.m)
        
class Ryd_pair (object):
    def __init__(self, atom1, atom2):
        """
        Define a pair of atom 1 & 2. Return E12, E12_Zeeman
        """
        if isinstance(atom1, Ryd_atom) & isinstance(atom2, Ryd_atom):
            self.atom1 = atom1
            self.atom2 = atom2
            self. E = atom1.En + atom2.En
            self.E_Zeeman = self.E + atom1.Zeemanshift + atom2.Zeemanshift
        else:
            raise Exception('input should be of Ryd_atom class')
            pass
            
    def __repr__(self):
        return '{0} and {1}'.format(self.atom1, self.atom2)
    def __eq__(self, other):
        return (self.atom1 == other.atom1) & (self.atom2 == other.atom2)
def pair_invert(pair):
    """ pair_invert(pair)
    Invert pair of atoms"""
    return Ryd_pair(pair.atom2, pair.atom1)
