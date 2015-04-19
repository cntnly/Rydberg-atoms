# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from matplotlib import *
from pylab import *
import builtins
#import numexpr as ne
import sys
path = 'C:/Users/r14/Documents/GitHub/test/python'
if path not in sys.path:
    sys.path.append(path)

from imp import reload

if __name__ == '__main__':
    import para
    from para import *
        
    reload(para)
    from para import *
#import constant
#from constant import * # define constant to be used

#reload(constant)
#from constant import *
# import radinte
import pyximport; pyximport.install()
try:
    from radinte_sqrt import *
except ImportError:
    raise Exception('No radinte module found')
try:
    from functions import *
except ImportError:
    raise Exception('No functions module found')
    


# Calculate Energy
atom_1 = Ryd_atom(n1, l1, m1)
atom_2 = Ryd_atom(n2, l2, m2)
pair_12 = Ryd_pair(atom_1, atom_2)

coef = a_0*a_0*e*e/(4*np.pi*epsilon_0*h)
test_term = coef/(R_test**3)
test_term = test_term*test_term

print(pair_12)
print('theta = {0} deg'.format(theta*180/pi))
print('B_field = {0} G'.format(Bfield*1e4))
print('F_field = {0} V/cm'.format(Ffield))
print('theta_F = {0} deg'.format(theta_F*180/pi))
print('phi_F = {0} deg'.format(phi_F*180/pi))

builtins.Ffield = Ffield*100 # V/m
coef_F = Ffield*e*a_0/h

from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure('config')
clf()
ax = fig.gca(projection='3d')
ax.plot([0,0],[0,0],[0,1],'red')
ax.plot([0,np.sin(theta)],[0,0],[0,np.cos(theta)],'-o')
ax.plot([0,Ffield*np.sin(theta_F)*np.cos(phi_F)],[0,Ffield*np.sin(theta_F)*np.sin(phi_F)],[0,Ffield*np.cos(theta_F)])
plt.show()

#check degeneracy
N_list = [pair_12]
#N_deg=[]
#for lA in np.arange(max(l1-1, 0), l1+1.1,1):
#    if lA < n1:
#        for mA in np.arange(m1-1, m1+ 1.1,1):
#            if np.abs(mA) <= lA:                 
#                try:
#                    atomA_temp = Ryd_atom(n1, lA, mA)
#                except Exception:
#                    print(lA, mA)
#                if np.abs(atomA_temp.E_Zeeman - atom_1.E_Zeeman) < 1e-10:
#                    for lB in np.arange(max(l2-1, 0), l2+1.1):
#                        if lB < n2:
#                            for mB in np.arange(m2-1, m2+ 1.1):
#                                if np.abs(mB) <= lB:
#                                    try:
#                                        atomB_temp = Ryd_atom(n2, lB, mB)
#                                    except Exception:
#                                        print (lB, mB)
#                                    if np.abs(atomB_temp.E_Zeeman - atom_2.E_Zeeman) < 1e-10:
#                                        #print (atomB_temp)
#                                        pair_temp = Ryd_pair(atomA_temp, atomB_temp)
#                                        if pair_temp not in (N_list + N_deg):
#                                            N_deg.append(pair_temp)
#deg=len(N_deg)
#print('Degenerated pairs ={0}'.format(len(N_deg)+1))
##Search 1st order coupling terms
#N_list1 = Search_VdW(N_list, N_list + N_deg, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)
#N_list += N_list1 + N_deg
##N_list1 = Search_VdW(N_list, N_list, test_term, Choice, delta_n_max, l1, l2, 2, m1,m2, 3)
##N_list += N_list1
#
#print ('N_list = {0}'.format(len(N_list)))
#       
#
#for pair in N_list:
#    if pair_invert(pair) not in N_list:
#        N_list.append(pair_invert(pair))
#print ('1st order terms: {0}'.format(len(N_list)))
#
#
## Search for temp coupled by Stark effect
#N_list_Stark=[]
#Choice_F = (coef_F*coef_F)/Choice_F
#for elm in N_list:
#    N_list_Stark_temp = Search_Stark(elm, N_list + N_list_Stark, Choice_F, delta_n_max)
#    N_list_Stark += N_list_Stark_temp
#
#print('1st order Stark terms: {0}'.format(len(N_list_Stark)))
#
#count_doub(N_list, N_list_Stark)

#============= Create base =================
def create_base(pair, Not_list, delta_n, delta_l, delta_m, delta_E):
    atom_1 = pair.atom1
    atom_2 = pair.atom2
    N_list = []
    for nA in arange(atom_1.n -delta_n, atom_1.n +delta_n+0.1, 1):
        for lA in arange(max(0, atom_1.l -delta_l), min(atom_1.l +delta_l +0.1,nA),1):
            for mA in arange(-lA, lA+0.1,1):
                if abs(mA -atom_1.m) <=delta_m:
                    try:
                        atomA_temp = Ryd_atom(nA,lA,mA)
                    except:
                        print(nA, lA,mA)                            
                 #   if abs(atomA_temp.En - atom_1.En) < delta_E:
                    for nB in arange(atom_2.n -delta_n, atom_2.n +delta_n+0.1, 1):
                        for lB in arange(max(0, atom_2.l -delta_l), min(atom_2.l +delta_l+.1,nB),1):
                            for mB in arange(-lB, lB+0.1,1):
                                if abs(mB -atom_2.m) <=delta_m:
                                    try:
                                        atomB_temp = Ryd_atom(nB,lB,mB)
                                    except:
                                        print(nB, lB,mB)                            
                                 #   if abs(atomB_temp.En - atom_2.En) < delta_E:
                                    pairAB_temp = Ryd_pair(atomA_temp, atomB_temp) 
                                    if abs(pairAB_temp.E_Zeeman - pair_12.E_Zeeman) < delta_E:
                                        if  pairAB_temp not in Not_list:
                                            N_list.append(pairAB_temp)
    return N_list
N_list1 = create_base(pair_12, N_list, 5, 20,2, 100e9/2)                            
Union_list = N_list + N_list1
#=============================

#Union_list = N_list + N_list_Stark # sure that elements of N_list are not in N_list2
for pair in Union_list:
    if pair_invert(pair) not in Union_list:
        Union_list.append(pair_invert(pair))
        
        
print ('Matrix size: {0}'.format(len(Union_list)))
#try:
#    raw_input("Press enter to continue")
#except:
#    input("Press Enter to continue")
#print('wait...')           
#sys.stdout.flush()

# sort list as energy order
Union_list = sorted(Union_list, key = lambda energy: energy.E_Zeeman)
#print(len(Union_list))                    
figure(1)
clf()
subplot(2,2,1)
plot([(elm.E_Zeeman - pair_12.E_Zeeman)*1e-9 for elm in Union_list],'-o')
xlabel('List N°')
ylabel('Rel. energy (GHz)')
# Create interaction Matrix
length = len(Union_list)
# create mask
V_VdW = 0j*np.empty((length,length))
V_R = 0j*np.empty_like(V_VdW)
V_A = 0j*np.empty_like(V_VdW)
V_Stark1 = 0j*np.empty_like(V_VdW)
V_Stark2 = 0j*np.empty_like(V_VdW)
# mask to calculate only lower triangular matrix
mask1 = np.tril(np.ones_like(V_VdW))!=0

#vectorize necessary functions
rad_vec = np.vectorize(radinte)
A_vec = np.vectorize(A_Integral)
A_Stark_vec = np.vectorize(A_Stark)

V_row = np.asarray([(elm.atom1.n, elm.atom1.l, elm.atom1.m, elm.atom1.E_radinte, elm.atom2.n, elm.atom2.l, elm.atom2.m, elm.atom2.E_radinte) for elm in Union_list])
Pair1_n1, Pair2_n1 = meshgrid(V_row[:,0],V_row[:,0])
Pair1_l1, Pair2_l1 = meshgrid(V_row[:,1],V_row[:,1])
Pair1_m1, Pair2_m1 = meshgrid(V_row[:,2],V_row[:,2])
Pair1_Erad1, Pair2_Erad1 = meshgrid(V_row[:,3],V_row[:,3])
Pair1_n2, Pair2_n2 = meshgrid(V_row[:,4],V_row[:,4])
Pair1_l2, Pair2_l2 = meshgrid(V_row[:,5],V_row[:,5])
Pair1_m2, Pair2_m2 = meshgrid(V_row[:,6],V_row[:,6])
Pair1_Erad2, Pair2_Erad2 = meshgrid(V_row[:,7],V_row[:,7])

# Saving memory
Pair1_n1, Pair1_l1, Pair1_m1, Pair1_Erad1 = Pair1_n1[mask1], Pair1_l1[mask1], Pair1_m1[mask1], Pair1_Erad1[mask1]
Pair1_n2, Pair1_l2, Pair1_m2, Pair1_Erad2 = Pair1_n2[mask1], Pair1_l2[mask1], Pair1_m2[mask1], Pair1_Erad2[mask1]
Pair2_n1, Pair2_l1, Pair2_m1, Pair2_Erad1 = Pair2_n1[mask1], Pair2_l1[mask1], Pair2_m1[mask1], Pair2_Erad1[mask1]
Pair2_n2, Pair2_l2, Pair2_m2, Pair2_Erad2 = Pair2_n2[mask1], Pair2_l2[mask1], Pair2_m2[mask1], Pair2_Erad2[mask1]

V_A[mask1] = A_vec(Pair1_l1, Pair1_m1, Pair1_l2, Pair1_m2, Pair2_l1, Pair2_m1, Pair2_l2, Pair2_m2, theta)
# Choose only elements different from 0
mask2 = V_A[mask1] !=0
V_R[V_A!=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)*rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
V_VdW[V_A!=0] = V_R[V_A!=0]*V_A[V_A!=0]
V_VdW = V_VdW + V_VdW.T.conj() - np.diag(V_VdW.diagonal())
V_VdW = V_VdW*1e-9 # convert to GHz.m^3

figure(1)
subplot(2,2,2)
imshow(np.real(V_VdW),)
colorbar()

V_R = np.zeros_like(V_VdW)
V_A = np.zeros_like(V_VdW)
V_A[mask1] = A_Stark_vec(Pair1_l1, Pair1_m1, Pair2_l1, Pair2_m1, theta_F, phi_F)*(Pair1_n2 == Pair2_n2) *(Pair1_l2 == Pair2_l2)*(Pair1_m2 == Pair2_m2)
mask2 = V_A[mask1] !=0
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)
    V_Stark1[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark1 = V_Stark1 + V_Stark1.T.conj() - np.diag(V_Stark1.diagonal())
    V_Stark1 = V_Stark1*1e-9
subplot(2,2,3)
imshow(np.imag(V_Stark1))
colorbar()

V_R = np.zeros_like(V_VdW)
V_A = np.zeros_like(V_VdW)
V_A[mask1] = A_Stark_vec(Pair1_l2, Pair1_m2, Pair2_l2, Pair2_m2, theta_F, phi_F)*(Pair1_n1 == Pair2_n1) *(Pair1_l1 == Pair2_l1)*(Pair1_m1 == Pair2_m1)
mask2 = V_A[mask1] !=0
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
    V_Stark2[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark2 = V_Stark2 + V_Stark2.T.conj() - np.diag(V_Stark2.diagonal())
    V_Stark2 = V_Stark2*1e-9
subplot(2,2,4)
imshow(np.imag(V_Stark2))
colorbar()

# Zero-th Energy
EI = np.diag(np.asarray([elm.E_Zeeman for elm in Union_list]))*1e-9

R_max = 30 #in um
R_min = 1.
#R = R_max
R_num = 200 # step in logarithm
index1 = Union_list.index(pair_12)
if pair_invert(pair_12) != pair_12:
    index2 = Union_list.index(pair_invert(pair_12))
else:
    index2 = index1
#index2 = Union_list.index(pair_invert(pair_34))
#index3 = Union_list.index(pair_invert(pair_56))
#index4 = Union_list.index(pair_invert(pair_78))

out_egr = np.empty((R_num, length))
out_coef = np.empty((2,R_num))
out_vector = np.empty((R_num,length, length))
R = np.logspace(log10(R_min), log10(R_max), num = R_num)
for i,elm in enumerate(R):
    #out_egr[i] = np.linalg.eigvalsh(EI + V_VdW* coef/(elm**3)) - pair_12.E_Zeeman
    out_egr[i] , out_vector[i] = np.linalg.eigh(EI + 1e18*V_VdW* coef/(elm**3) + coef_F*(V_Stark1 + V_Stark2))
  #  print(out_egr[np.iscomplex(out_egr[i])])
    out_egr[i] = out_egr[i] - pair_12.E_Zeeman*1e-9
    out_coef[0,i] = np.argmax(abs(out_vector[i][index1,:]))
    out_coef[1,i] = np.argmax(abs(np.delete(out_vector[i][index2,:],out_coef[0,i])))
    if out_coef[1, i] >= out_coef[0,i]:
        out_coef[1,i] +=1

    
out_egr1 = [out_egr[i, out_coef[0,i]] for i in range(R_num)]
out_egr2 = [out_egr[i, out_coef[1,i]] for i in range(R_num)]



figure(2)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
loglog(R, asarray(out_egr), R, asarray(out_egr1), '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')
figure(5);clf()
loglog(R, -asarray(out_egr), R, -asarray(out_egr1), '+')
#semilogx(R, out_egr[:,index3])
#semilogx(R, out_egr[:,index4])
figure(3)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
semilogx(R, out_egr, R, out_egr1, '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')

#fit data
def VdW(r,C):
    return C*(r**-6)
def pow_fit(r,C,n):
    return C*(r**-6)

#popt,pcov = curve_fit(VdW, R[100:], abs(out_egr1[100:]))
#popt2,pcov = curve_fit(VdW, R[80:], abs(out_egr2[80:]))
#figure(4)
#clf()
#loglog(R, VdW(R, popt), R,abs(out_egr1))
#loglog(R, VdW(R, popt2),R,abs(out_egr2))

from scipy.optimize import curve_fit

offset, off_vec = np.linalg.eigh(EI + coef_F*(V_Stark1 + V_Stark2))
offset = offset - pair_12.E_Zeeman*1e-9
i = np.argmax(abs(off_vec[index1,:]))
offset1 = offset[i]
k = np.argmax(abs(np.delete(off_vec[index2,:], i)))
if k >= i:
    k +=1
offset2 = offset[k]
popt1,pcov1 = curve_fit(pow_fit, R[75:135], out_egr1[75:135]-offset1, p0=(100,6))
#popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:]-offset1, p0=(100,6))
figure(4)
clf()
loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1-offset1), '+')
#loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1-offset1), '+')
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')
print('C6 = {0} GHz.um^6'.format(popt1))
print('R_1 MHz = {0} um'.format(R[np.max(np.where(abs(out_egr1-offset1)>1e-3))]))

figure(6);clf()
#semilogx(R, out_vector[:, index1, out_coef[0,-1]]**2)
semilogx(R, asarray([out_vector[i, index1, out_coef[0,i]] for i in range(R_num)])**2)

#popt3,pcov3 = curve_fit(VdW, R[100:], abs(out_egr3[100:]- pair_56.E_Zeeman+ pair_12.E_Zeeman))
#popt4,pcov4 = curve_fit(VdW, R[100:], abs(out_egr4[100:]- pair_78.E_Zeeman+ pair_12.E_Zeeman))
if index2 != index1:
  #  out_egr2 = np.asarray([out_egr[i, out_coef[1,i]] for i in range(R_num)])
    #out_egr2 = out_egr[:, index2]
    figure(2)
    loglog(R, asarray(out_egr), R, asarray(out_egr2), '+')
    figure(5)
    loglog(R, -asarray(out_egr), R, -asarray(out_egr2), '+')
    figure(3)
    semilogx(R, out_egr2, 'wo')
    popt2,pcov2 = curve_fit(pow_fit, R[100:], abs(out_egr2[100:]-offset2), p0=(100,6))
    figure(4)
    loglog(R, abs(pow_fit(R, *popt2)), R,abs(out_egr2-offset2),'wo')
    print(popt2)
  





#popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:])
#figure(5)
#clf()
#loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1))
#print(popt1)