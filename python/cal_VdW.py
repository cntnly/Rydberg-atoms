# -*- coding: utf-8 -*-
from __future__ import division
import numpy as np
from matplotlib import *
from pylab import *
import builtins
#import numexpr as ne
import sys, os
path = os.path.dirname(os.path.abspath('__file__'))
path= '/home/2kome/Desktop/testdeck/test1/python'
#path = 'C:/Users/r14/Documents/GitHub/test/python'
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
print('B_field effective = {0} G'.format(Bfield*1e4))
print('F_field = {0} V/cm'.format(Ffield))
print('theta_B = {0} deg'.format(theta_B*180/pi))
print('phi_B = {0} deg'.format(phi_B*180/pi))

builtins.Ffield = Ffield*100 # V/m
coef_F = Ffield*e*a_0/h

Bz = Bfield*np.cos(theta_B)*mu_B/h
Bx = Bfield*np.sin(theta_B)*np.cos(phi_B)*mu_B/h
By = Bfield*np.sin(theta_B)*np.sin(phi_B)*mu_B/h

from mpl_toolkits.mplot3d import Axes3D
figure(0)
clf()
ax=subplot(2,2,1, projection ='3d')
#ax = fig.gca(projection='3d')
plot([0,0],[0,0],[0,1.2],'red') # F_field
plot([0,np.sin(theta)],[0,0],[0,np.cos(theta)],'-o') # atoms
plot([0,0+Bx],[0,0+By],[0,0+Bz]) # B_field
ax.set_xticks([]);ax.set_yticks([]);ax.set_zticks([]);

#check degeneracy
N_list = [pair_12]


#============= Create base =================

N_list1 = create_base(pair_12, N_list,2, 20,2, 100e9/2)                            
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
Union_list = sorted(Union_list, key = lambda energy: energy.E)
#print(len(Union_list))  
if __name__ == '__main__':             
    figure(1)
    clf()
    subplot(2,2,1)
    plot([(elm.E - pair_12.E)*1e-9 for elm in Union_list],'-o')
    xlabel('List N°')
    ylabel('Rel. energy (GHz)')
# Create interaction Matrix
length = len(Union_list)

# create mask
V_VdW = np.zeros((length,length))
V_R = np.zeros_like(V_VdW)
V_A = np.zeros_like(V_VdW)
V_Stark1 = np.zeros_like(V_VdW)
V_Stark2 = np.zeros_like(V_VdW)
V_Zeeman = np.zeros_like(V_VdW)

# mask to calculate only lower triangular matrix
mask1 = (np.tril(np.ones_like(V_VdW))!=0)

#vectorize necessary functions
rad_vec = np.vectorize(radinte)
A_vec = np.vectorize(A_Integral)
A_Stark_vec = np.vectorize(A_Stark)
A_Zeeman_vec = np.vectorize(A_Zeeman)

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
mask2 = (V_A[mask1] !=0)
V_R[V_A!=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)*rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
V_VdW[V_A!=0] = V_R[V_A!=0]*V_A[V_A!=0]
V_VdW = V_VdW + V_VdW.T - np.diag(V_VdW.diagonal())
V_VdW = V_VdW*1e-9 # convert to GHz.m^3

if __name__ == '__main__':             
    figure(1)
    subplot(2,2,2)
    imshow(V_VdW)
    colorbar()

V_R = np.zeros_like(V_VdW)
V_A = np.zeros_like(V_VdW)
V_A[mask1] = A_Stark_vec(Pair1_l1, Pair1_m1, Pair2_l1, Pair2_m1)*(Pair1_n2 == Pair2_n2) *(Pair1_l2 == Pair2_l2)*(Pair1_m2 == Pair2_m2)
mask2 = (V_A[mask1] !=0)
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad1[mask2], Pair1_l1[mask2], Pair2_Erad1[mask2], Pair2_l1[mask2], 1)
    V_Stark1[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark1 = V_Stark1 + V_Stark1.T - np.diag(V_Stark1.diagonal())
    V_Stark1 = V_Stark1*1e-9
if __name__ == '__main__':             
    subplot(2,2,3)
    imshow(V_Stark1)
    colorbar()

V_R = np.zeros_like(V_VdW)
V_A = np.zeros_like(V_VdW)
V_A[mask1] = A_Stark_vec(Pair1_l2, Pair1_m2, Pair2_l2, Pair2_m2)*(Pair1_n1 == Pair2_n1) *(Pair1_l1 == Pair2_l1)*(Pair1_m1 == Pair2_m1)
mask2 = (V_A[mask1] !=0)
if len(V_A[V_A!=0])!=0:
    V_R[V_A !=0] = rad_vec(Pair1_Erad2[mask2], Pair1_l2[mask2], Pair2_Erad2[mask2], Pair2_l2[mask2], 1)
    V_Stark2[V_A !=0] = V_R[V_A!=0]*V_A[V_A!=0]
    V_Stark2 = V_Stark2 + V_Stark2.T - np.diag(V_Stark2.diagonal())
    V_Stark2 = V_Stark2*1e-9
if __name__ == '__main__':             
    subplot(2,2,4)
    imshow(V_Stark2)
    colorbar()

# Zeeman
V_A = 0j*np.zeros_like(V_VdW)
V_A[mask1] = A_Zeeman_vec(Pair1_l1, Pair1_m1, Pair2_l1, Pair2_m1, Bx, By, Bz)*(Pair1_n1 == Pair2_n1)*(Pair1_n2 == Pair2_n2) *(Pair1_l2 == Pair2_l2)*(Pair1_m2 == Pair2_m2)
V_Zeeman += V_A

V_A = 0j*np.zeros_like(V_VdW)
V_A[mask1] = A_Zeeman_vec(Pair1_l2, Pair1_m2, Pair2_l2, Pair2_m2, Bx, By, Bz)*(Pair1_n2 == Pair2_n2) *(Pair1_n1 == Pair2_n1) *(Pair1_l1 == Pair2_l1)*(Pair1_m1 == Pair2_m1)
V_Zeeman += V_A

V_Zeeman = V_Zeeman + V_Zeeman.T - np.diag(V_Zeeman.diagonal())
V_Zeeman = V_Zeeman*1e-9
    
# Zero-th Energy
EI = np.diag(np.asarray([elm.E for elm in Union_list]))*1e-9

R_max = 200 #in um
R_min = 1.
#R = R_max
R_num = 200 # step in logarithm
index1 = Union_list.index(pair_12)
if pair_invert(pair_12) != pair_12:
    index2 = Union_list.index(pair_invert(pair_12))
else:
    index2 = index1

out_egr = np.empty((R_num, length))
out_coef = np.empty((2,R_num))
out_vector = np.empty((R_num,length, length))
R = np.logspace(log10(R_min), log10(R_max), num = R_num)
for i,elm in enumerate(R):
    #out_egr[i] = np.linalg.eigvalsh(EI + V_VdW* coef/(elm**3)) - pair_12.E_Zeeman
    out_egr[i] , out_vector[i] = np.linalg.eigh(EI + 1e18*V_VdW* coef/(elm**3) + coef_F*(V_Stark1 + V_Stark2)+ V_Zeeman)
  #  print(out_egr[np.iscomplex(out_egr[i])])
    out_egr[i] = out_egr[i] - pair_12.E*1e-9
    
    if i == 0:
        out_coef[0,i] = np.argmax(abs(out_vector[i][index1,:]))
    else:
        out_coef[0,i] = np.argmax(abs(out_vector[i][index1,:]))
        if (out_coef[0,i]!= out_coef[0,i-1]) and abs(out_vector[i][index1, out_coef[0,i-1]]**2 -out_vector[i][index1, out_coef[0,i]]**2) <0.001:
            #print(i, out_coef[0,i-1], out_vector[i][index1, out_coef[0,i-1]]**2, out_coef[0,i], out_vector[i][index1, out_coef[0,i]]**2 )
            out_coef[0,i] = out_coef[0,i-1]
    out_coef[1,i] = np.argmax(abs(np.delete(out_vector[i][index2,:],out_coef[0,i])))
    if out_coef[1, i] >= out_coef[0,i]:
        out_coef[1,i] +=1

    
out_egr1 = [out_egr[i, out_coef[0,i]] for i in range(R_num)]
out_egr2 = [out_egr[i, out_coef[1,i]] for i in range(R_num)]

out_vec1 = np.asarray([out_vector[i,:,out_coef[0,i]] for i in range(R_num)])
if index1!=index2:
    out_vec2 = np.asarray([out_vector[i,:,out_coef[1,i]] for i in range(R_num)])


#figure(2)
#clf()
##semilogx(R, out_egr1)
##semilogx(R, out_egr2)
#loglog(R, asarray(out_egr), R, asarray(out_egr1), '+')
#xlabel('$R (\mu$m)')
#ylabel('Rel. energy (GHz)')

figure(0);subplot(2,2,2);
pcolor(R,np.arange(length+1)-0.5,out_vec1.T**2)
xlim(R_min, R_max), ylim(ymin = max(index1 -100.5,-0.5),ymax= min(index1+100.5,length-.5))
xscale('log')
xticks([])
tick_params(which='both', direction='out')

#loglog(R, -asarray(out_egr), R, -asarray(out_egr1), '+')
#semilogx(R, out_egr[:,index3])
#semilogx(R, out_egr[:,index4])
if __name__ == '__main__':
    figure(2)
    clf()
    #semilogx(R, out_egr1)
    #semilogx(R, out_egr2)
    semilogx(R, out_egr, R, out_egr1, '+')
    xlim(R_min, R_max)
    xlabel('$R (\mu$m)')
    ylabel('Rel. energy (GHz)')

figure(3)
clf()
#semilogx(R, out_egr1)
#semilogx(R, out_egr2)
semilogx(R, out_egr, R, out_egr1, '+')
xlim(R_min, R_max)
ylim(out_egr1[-1]-0.05, out_egr1[-1]+0.1)
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')

#fit data
def VdW(r,C):
    return C*(r**-6)
def pow_fit(r,C,p):
    return C*(r**-p)

#popt,pcov = curve_fit(VdW, R[100:], abs(out_egr1[100:]))
#popt2,pcov = curve_fit(VdW, R[80:], abs(out_egr2[80:]))
#figure(4)
#clf()
#loglog(R, VdW(R, popt), R,abs(out_egr1))
#loglog(R, VdW(R, popt2),R,abs(out_egr2))

from scipy.optimize import curve_fit

offset, off_vec = np.linalg.eigh(EI + coef_F*(V_Stark1 + V_Stark2)+V_Zeeman)
offset = offset - pair_12.E*1e-9
i = np.argmax(abs(off_vec[index1,:]))
offset1 = offset[i]
k = np.argmax(abs(np.delete(off_vec[index2,:], i)))
if k >= i:
    k +=1
offset2 = offset[k]
popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:]-offset1, p0=(100,6))
#popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:]-offset1, p0=(100,6))
p1 = round(popt1[1])
popt1,pcov1 = curve_fit(lambda R, C: pow_fit(R,C,p1), R[100:], out_egr1[100:]-offset1, p0=popt1[0])

figure(0);
subplot(2,2,3)
loglog(R,abs(out_egr1-offset1), '+',R, abs(pow_fit(R, popt1, p1)))
#a,b =ax2.get_ylim()
loglog(R, R**-3, 'g--', R, 100*R**-6,'b--',linewidth=0.5)
#loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1-offset1), '+')
xlim(R_min, R_max)
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')
#print('R_1 MHz = {0} um'.format(R[np.max(np.where(abs(out_egr1-offset1)>1e-3))]))

figure(0);subplot(2,2,4)
#semilogx(R, out_vector[:, index1, out_coef[0,-1]]**2)
semilogx(R, asarray([out_vector[i, index1, out_coef[0,i]] for i in range(R_num)])**2)
xlim(R_min, R_max)
xlabel('$R (\mu$m)')
ylabel('a.u.')
#popt3,pcov3 = curve_fit(VdW, R[100:], abs(out_egr3[100:]- pair_56.E_Zeeman+ pair_12.E_Zeeman))
#popt4,pcov4 = curve_fit(VdW, R[100:], abs(out_egr4[100:]- pair_78.E_Zeeman+ pair_12.E_Zeeman))
if index2 != index1:
  #  out_egr2 = np.asarray([out_egr[i, out_coef[1,i]] for i in range(R_num)])
    #out_egr2 = out_egr[:, index2]
#    figure(2)
#    loglog(R, asarray(out_egr), R, asarray(out_egr2), '+')
#    figure(5)
#    loglog(R, -asarray(out_egr), R, -asarray(out_egr2), '+')
#    figure(3)
#    semilogx(R, out_egr2, 'wo')
#    popt2,pcov2 = curve_fit(pow_fit, R[100:], abs(out_egr2[100:]-offset2), p0=(100,6))
#    figure(4)
#    loglog(R, abs(pow_fit(R, *popt2)), R,abs(out_egr2-offset2),'wo')
#    print(popt2)
    popt2,pcov2 = curve_fit(pow_fit, R[100:], out_egr2[100:]-offset2, p0=(1,3))
    #popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:]-offset1, p0=(100,6))
    p2 = round(popt2[1])
    popt2,pcov2 = curve_fit(lambda R, C: pow_fit(R,C,p2), R[100:], out_egr2[100:]-offset2, p0=popt2[0])
    
    figure(0)
    subplot(2,2,3)
    loglog(R,abs(out_egr2-offset2), 'wo',R, abs(pow_fit(R, popt2, p2)))
#    print('C{0:0.0f} = {1} GHz.um^{0:0.0f}'.format(p2,popt2))
    figure(3)
    semilogx(R, out_egr2, 'wo')  
    if __name__ == '__main__':
        figure(2)
        semilogx(R, out_egr2, 'wo')
    
#    for i in range(R_num):
#        x=where(out_vec2[i]**2>0.05)
#        if len(x[0]) >1:
#            print('At R = {0} um, prop > 5% are:'.format(R[i]))
#            for elm in x[0]:
#                print('pair of {0}, at index{1}'.format(Union_list[elm],elm))
#     #       print([(elm,Union_list[elm]) for elm in x[0]])
#            break

show()
print('C{0:0.0f} = {1} GHz.um^{0:0.0f}'.format(p1,popt1))
if index1!=index2:
    print('C{0:0.0f} = {1} GHz.um^{0:0.0f}'.format(p2,popt2))
  
for i in range(R_num-1, -1,-1):
    count = 0
    x=where(out_vec1[i]**2>0.05)
    if index1 in x[0]:
        count +=1
    if (index2!= index1) &(index2 in x[0]):
        count +=1
    if len(x[0]) >count:
        print('At R = {0} um, prop > 5% are:'.format(R[i]))
        for elm in x[0]:
            print('pair of {0}, at index {1}'.format(Union_list[elm],elm))
 #       print([(elm,Union_list[elm]) for elm in x[0]])
        break


#popt1,pcov1 = curve_fit(pow_fit, R[100:], out_egr1[100:])
#figure(5)
#clf()
#loglog(R, abs(pow_fit(R, *popt1)), R,abs(out_egr1))
#print(popt1)
def index(arg):
    return max(where(R<=arg)[0])
def doubled(r,C3,C6):
    return C3*(r**-3) + C6*(r**-6)    
popt1,pcov1 = curve_fit(doubled, R[100:], out_egr1[100:]-offset1, p0=(-1,1))
figure(6); clf();
loglog(R,abs(out_egr1-offset1), '+',R, abs(doubled(R, *popt1)))
loglog(R, R**-3, 'g--', R, 100*R**-6,'b--',linewidth=0.5)
xlim(R_min, R_max)
xlabel('$R (\mu$m)')
ylabel('Rel. energy (GHz)')

print(popt1)
