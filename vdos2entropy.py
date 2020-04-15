#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.font_manager
import scipy
import math as m
import sys, os

#############################
# Parse cmd line arguments: #
#############################

if (len(sys.argv) - 1) != 7:
    print("Illegal number of command line arguments!")
    print("Required files:")
    print("  (1) working directory")
    print("  (2) file name COM VDOS")
    print("  (3) file name ROT VDOS (x component)")
    print("  (4) file name ROT VDOS (y component)")
    print("  (5) file name ROT VDOS (z component)")
    print("  (6) width for moving averages used in smoothing")
    print("  (7) base name of output files")
    sys.exit(2)

#Use LaTeX fonts
rc('text.latex', preamble=r'\usepackage{cmbright}')
mpl.rcParams['text.usetex'] = True

# Physical constants
kB = 1.3806488*10**(-23)          #Boltzmann's constant
T = 300                           #Temperature (K)
nDOF = 3.0                        #Degrees of freedom
h = 6.62606957*10**(-34)          #Planck's constant
Na = 6.022*10**23                 #Avogadro number
mO = 15.9994                      #Oxygen mass (amu)
mH = 1.0079                       #Hydrogen mass (amu)
mH2O = (mO + 2.0*mH)*0.001/Na     #Water mass (amu)
OH = 0.9572                       #O-H bond distance (Angström)
HOH = 104.52                      #H-O-H angle (°)

# Reference values for bulk water at T=300 K using the 2PT method.
# (taken from Fisette, Päslack et al., JACS (2016), DOI: 10.1021/jacs.6b07005)
ref_trans = 46.297
ref_rot = 9.303
ref_tot = ref_trans + ref_rot
ref_exp = 69.95

##############
# Functions: #
##############

def check_file(file_name):
    try:
        with open(file_name) as f:
            f.readlines()
            exit_status = 0
    except FileNotFoundError:
        print("file " + file_name + " does not exist!")
        exit_status = 1
        #sys.exit()
    return exit_status

def smooth_spectra(spec,nBins):
    freq = spec[0]

    dim = int(len(freq)/nBins)
    coarseSpec = np.zeros([2,dim])

    n=0
    for i in range(0,len(freq),nBins):
        for j in range(i,i+nBins,1):
            coarseSpec[0][n] += spec[0][j]
            coarseSpec[1][n] += spec[1][j]
        coarseSpec[0][n] /= nBins
        coarseSpec[1][n] /= nBins
        n += 1

    return coarseSpec

def COMEntropySolidQM(freq,specdens,temp):

    factor = h * (freq / 33.3564) * (10**12) / (kB*temp)
    s = specdens * (factor / (np.exp(factor) - 1.0) - np.log(1.0 - np.exp(-factor)))
    s = kB * s * (df / 33.3564) * (10**12)
    return s

def ROTEntropySolidQM(freq,specdens,temp):

    factor = h * (freq / 33.3564) * (10**12) / (kB*temp)
    s = specdens * (factor / (np.exp(factor) - 1.0) - np.log(1.0 - np.exp(-factor)))
    s = kB * s * (df / 33.3564) * (10**12)
    return s

#########################
# Read the input files: #
#########################

dir = sys.argv[1]
fnCOM = sys.argv[2]
fnROT = [sys.argv[3], sys.argv[4], sys.argv[5]]

nSim = 1
for subdir in os.listdir(dir):
    if "nve-" in subdir:
        nSim = nSim + 1

print(dir)
print("nSim = ",str(nSim-1))

n=0
for subdir in os.listdir(dir):
    exit_status = 0
    if "nve-" in subdir:
        current_dir = dir + "/" + subdir
        fn1 = current_dir + "/" + fnCOM
        fn2 = [ current_dir + "/" + fnROT[0],
                current_dir + "/" + fnROT[1],
                current_dir + "/" + fnROT[2]]

        #check if files exist, if not continue to next directory
        exit_status += check_file(fn1)
        for i in range(len(fn2)):
            exit_status += check_file(fn2[i])
        if exit_status != 0:
            continue

        x1,y1 = np.genfromtxt(fn1, unpack=True)
        x2,y2 = np.empty([3,len(x1)]), np.empty([3,len(x1)])
        for i in range(3):
            x2[i], y2[i] = np.genfromtxt(fn2[i], unpack=True)
        if n==0:
            vdosCOM = np.zeros([2,len(x1)])
            vdosROT = np.zeros([3,2,len(x1)])
        for j in range(len(x1)):
            vdosCOM[0][j] = x1[j]
            vdosCOM[1][j] += y1[j]
            for i in range(3):
                vdosROT[i][0][j] = x2[i,j]
                vdosROT[i][1][j] += y2[i,j]
        n = n + 1

# Averaging over simulations and smoothing spectra by moving averages
vdosCOM[1] = vdosCOM[1] / n
nbins = int(sys.argv[6])
vdos = smooth_spectra(vdosCOM,nbins)

rdos = np.zeros([3,2,round(len(x1)/nbins)])
for i in range(3):
    vdosROT[i][1] = vdosROT[i][1] / n
    rdos[i] = smooth_spectra(vdosROT[i],nbins)

df = vdos[0,1] - vdos[0,0]
print("frequency spacing: " + str(df))

##########################
# Translational entropy: #
##########################

# normalize vdos such that integral gives 3N-6 DOF
integralCOM = np.sum(vdos[1]) * (df / 33.3564 * 10**12)
vdosCOMnorm = vdos[1] / integralCOM * nDOF

Scom = np.zeros([len(vdos[1])])
for i in range(1,len(vdos[1])):
    if i==1:
        Scom[i] = COMEntropySolidQM(vdos[0,i],vdosCOMnorm[i],T) * Na
    else:
        Scom[i] = Scom[i-1] + COMEntropySolidQM(vdos[0,i],vdosCOMnorm[i],T) * Na

error_trans = (1.0-Scom[-1]/ref_trans)*100

print("entropy (trans): " + str(Scom[-1]))
print("  error vs 2PT: " +str(error_trans)+"%")

#######################
# Rotational entropy: #
#######################

#Define coordinate system based on internal water coordinates
waterCoord = [np.zeros([3]),
              np.dot([
              	      [np.cos(np.deg2rad(0.5 * HOH)),0,np.sin(np.deg2rad(0.5 * HOH))],
              	      [0,1,0],
              	      [-np.sin(np.deg2rad(0.5 * HOH)),0,np.cos(np.deg2rad(0.5 * HOH))]
              	     ],
                     [0,0,-OH]),
              np.dot([
              	      [np.cos(-np.deg2rad(0.5 * HOH)),0,np.sin(-np.deg2rad(0.5 * HOH))],
              	      [0,1,0],
              	      [-np.sin(-np.deg2rad(0.5 * HOH)),0,np.cos(-np.deg2rad(0.5 * HOH))]
              	      ],
                     [0,0,-OH])]

#Shift w.r.t to center of mass
masses = np.zeros([3])
masses = mO,mH,mH
waterCOM = np.matmul(masses,waterCoord)/(mO + 2*mH)

waterCoord = waterCoord - waterCOM

#Define moments of inertia around molecular axes
IA = (  mO * np.linalg.norm([waterCoord[0, 1], waterCoord[0, 2]])**2
      + mH * np.linalg.norm([waterCoord[1, 1], waterCoord[1, 2]])**2
      + mH * np.linalg.norm([waterCoord[2, 1], waterCoord[2, 2]])**2)/(Na*1000*10**20)

IB = (  mO * np.linalg.norm([waterCoord[0, 0], waterCoord[0, 2]])**2
      + mH * np.linalg.norm([waterCoord[1, 0], waterCoord[1, 2]])**2
      + mH * np.linalg.norm([waterCoord[2, 0], waterCoord[2, 2]])**2)/(Na*1000*10**20)

IC = (  mO * np.linalg.norm([waterCoord[0, 0], waterCoord[0, 1]])**2
      + mH * np.linalg.norm([waterCoord[1, 0], waterCoord[1, 1]])**2
      + mH * np.linalg.norm([waterCoord[2, 0], waterCoord[2, 1]])**2)/(Na*1000*10**20)

inertia = [IA, IB, IC]

#Weight vdos by moments of intertia
vdosROTweighted = np.empty([3,len(rdos[1,1])])
for i in [0,1,2]:
    vdosROTweighted[i] = rdos[i,1] / inertia[i]

#Average rotational vdos over all components
vdosROTmean = np.mean(vdosROTweighted, axis=0)

#normalize vdos such that integral gives 3N-6 DOF
integralROT = np.sum(vdosROTmean) * (df / 33.3564 * 10**12)
vdosROTnorm = vdosROTmean / integralROT * nDOF

Srot = np.zeros([len(vdosROTnorm)])
for i in range(1,len(vdosROTnorm)):
    if i==1:
        Srot[i] = ROTEntropySolidQM(vdos[0,i],vdosROTnorm[i],T) * Na
    else:
        Srot[i] = Srot[i-1] + ROTEntropySolidQM(vdos[0,i],vdosROTnorm[i],T) * Na

error_rot = (1.0-Srot[-1]/ref_rot)*100

print("entropy (rot): " + str(Srot[-1]))
print("  error vs 2PT: " +str(error_rot)+"%")

#Compute total entropy
print("total entropy: " + str(int(Srot[-1] + Scom[-1])) + " J/(mol K)")

#######################
# Save data to files: #
#######################

base_name = sys.argv[7]

outfile = open(base_name + "_entropy.dat","w+")
for i in range(len(Scom)):
     outfile.write("%f\t%f\t%f\t%f\n" % (vdos[0,i], Scom[i], Srot[i], Scom[i]+Srot[i]))
outfile.close()

outfile = open(base_name + "_vdos.dat","w+")
for i in range(len(vdos[0])):
     outfile.write("%f\t%f\t%f\t%f\n" % (vdos[0,i], vdosCOMnorm[i], vdosROTnorm[i], vdosCOMnorm[i]+vdosROTnorm[i]))
outfile.close()

###################
# Plot entropies: #
###################

error_tot = (1.0 - (Srot[-1] + Scom[-1])/ref_tot) * 100.0

fig, ax = plt.subplots(figsize=(4, 3))
x1, y1 = [0, 1200], [ref_tot, ref_tot]
x2, y2 = [0, 1200], [ref_rot , ref_rot]
x3, y3 = [0, 1200], [ref_trans, ref_trans]
ax.plot(vdos[0], Srot,color='r',label='SRE rot')
ax.plot(vdos[0], Scom,color='b',label='SRE trans')
ax.plot(vdos[0], Scom+Srot,color='gray',label='SRE tot')
ax.plot(x1,y1,color='gray',linestyle='--',label='total bulk (2PT)')
ax.plot(x2,y2,color='r',linestyle='--',label='rot bulk (2PT)')
ax.plot(x3,y3,color='b',linestyle='--',label='trans bulk (2PT)')
ax.set(xlabel=r'Frequency ($\mathrm{cm^{-1}}$)', ylabel=r'Entropy ($\mathrm{J/K\cdot mol}$)')
ax.text(100,12,'deviation from bulk: '+ str(int(ref_tot - (Srot[-1] + Scom[-1]))) + r'$\mathrm{J/K\cdot mol}$')
plt.legend(loc='center right',fontsize=7)
plt.xlim(0, 1200)
plt.ylim(0, 60)
plt.tight_layout()
plt.savefig(base_name + "_entropy_SRE" + ".pdf")
#plt.show()

#########################
# Plot normalized VDOS: #
#########################

fig, ax = plt.subplots(figsize=(4, 3))
ax.plot(vdos[0],(vdosROTnorm+vdosCOMnorm)*10**13, color='black', label='total')
ax.plot(vdos[0], vdosROTnorm*10**13, color='red', label='rotation')
ax.plot(vdos[0], vdosCOMnorm*10**13, color='blue', label='translation')
ax.set(xlabel=r'Frequency ($\mathrm{cm^{-1}}$)', ylabel=r'Spectral density (arb. u.)')
plt.legend(loc='upper right',fontsize=9)
plt.xlim(0, 1200)
plt.ylim(0, 10)
plt.tight_layout()
plt.savefig(base_name + "_vdos" + ".pdf")
#plt.show()
