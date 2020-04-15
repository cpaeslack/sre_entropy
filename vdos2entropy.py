#!/usr/bin/env python

import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib.font_manager
import math as m
import sys, os
import smooth_spectrum as s

##############
# Functions: #
##############

def check_file(file_name):
  try:
    f = open("filename.txt")
  except FileNotFoundError:
    print(file_name + " not accessible!")
  finally:
    f.close()

def entropySolidQM(freq,specdens,temp):

    factor = h * (freq / 33.3564) * (10**12) / (kB*temp)
    s = specdens * (factor / (np.exp(factor) - 1.0) - np.log(1.0 - np.exp(-factor)))
    s = kB * s * (df / 33.3564) * (10**12)
    return s

def print_arguments():
  print("Required arguments:")
  print("  (1) working directory")
  print("  (2) file name COM VDOS")
  print("  (3) file name ROT VDOS (1st component)")
  print("  (4) file name ROT VDOS (2nd component)")
  print("  (5) file name ROT VDOS (3rd component)")
  print("  (6) width to be used for smoothing")
  print("  (7) base name of output files")

#############################
# Parse cmd line arguments: #
#############################

if sys.argv[1] == "--help":
  print("This script computes the configurational entropy")
  print("of water based on the spectral density,")
  print("as previously described by Schlitter & Massarczyk")
  print("(http://arxiv.org/abs/1909.04726).")
  print("")
  print("Input files need to be spectral densities computed")
  print("from the velocity auto-correlation function, i.e.")
  print("  - translational DoS (from mass-weighted cartesian velocities)")
  print("  - rotational DoS (individual 3 angular velocities along")
  print("    the molecular axes of water molecules")
  print("")
  print_arguments()
  sys.exit(0)

if (len(sys.argv) - 1) != 7:
  print("Illegal number of arguments!")
  print_arguments()
  sys.exit(0)

###########################
# Settings and constants: #
###########################

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

#########################
# Read the input files: #
#########################

dir = sys.argv[1]
fnCOM = sys.argv[2]
fnROT = [sys.argv[3], sys.argv[4], sys.argv[5]]

#check if files exist, if not continue to next directory
check_file(fn1)
for f in fn2:
  check_file(f)


x, y = np.genfromtxt(fn1, unpack=True)

vdosCOM = np.zeros([2,len(x)])
vdosROT = np.zeros([3,2,len(x)])

vdosCOM[0] = np.copy(x)
vdosCOM[1] = np.copy(y)

for i in range(3):
  x, y = np.genfromtxt(fn2[i], unpack=True)
  vdosROT[i][0] = np.copy(x[i])
  vdosROT[i][1] = np.copy(y[i])

nfreq = len(vdosCOM[1])
print("nfreq = " + str(nfreq))

# Averaging over simulations and gaussian smoothing
win_size = int(sys.argv[6])
print("win_size = " + str(win_size))

vdos = np.zeros([2,nfreq])
vdos[0] = np.copy(vdosCOM[0])
vdos[1] = np.copy(s.smooth_gaussian(vdosCOM[1],win_size))

rdos = np.zeros([3,2,nfreq])
for i in range(3):
  vdosROT[i][1] = vdosROT[i][1]
  rdos[i][0] = np.copy(vdosROT[i][0])
  rdos[i][1] = np.copy(s.smooth_gaussian(vdosROT[i][1],win_size))

df = vdos[0,1] - vdos[0,0]
print("df: " + str(df))

##########################
# Translational entropy: #
##########################

# normalize VDoS such that integral gives 3N-6 DOF
integralCOM = np.sum(vdos[1]) * (df / 33.3564 * 10**12)
vdosCOMnorm = vdos[1] / integralCOM * nDOF

Scom = np.zeros([len(vdos[1])])
for i in range(1,len(vdos[1])):
  if i==1:
    Scom[i] = entropySolidQM(vdos[0,i],vdosCOMnorm[i],T) * Na
  else:
    Scom[i] = Scom[i-1] + entropySolidQM(vdos[0,i],vdosCOMnorm[i],T) * Na

error_trans = (1.0-Scom[-1]/ref_trans)*100

print("entropy (trans): " + str(Scom[-1]))

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

#Shift w.r.t. to center of mass
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

#Weigh VDoS by moments of intertia
vdosROTweighted = np.empty([3,len(rdos[1,1])])
for i in [0,1,2]:
  vdosROTweighted[i] = rdos[i,1] / inertia[i]

#Average rotational VDoS over all components
vdosROTmean = np.mean(vdosROTweighted, axis=0)

#normalize VDoS such that integral gives 3N-6 DOF
integralROT = np.sum(vdosROTmean) * (df / 33.3564 * 10**12)
vdosROTnorm = vdosROTmean / integralROT * nDOF

Srot = np.zeros([len(vdosROTnorm)])
for i in range(1,len(vdosROTnorm)):
  if i==1:
    Srot[i] = entropySolidQM(vdos[0,i],vdosROTnorm[i],T) * Na
  else:
    Srot[i] = Srot[i-1] + entropySolidQM(vdos[0,i],vdosROTnorm[i],T) * Na

error_rot = (1.0-Srot[-1]/ref_rot)*100

print("entropy (rot): " + str(Srot[-1]))

#Compute total entropy
print("total entropy: " + str(int(Srot[-1] + Scom[-1])) + " J/(mol K)")

############################
# Compute entropy density  #
############################

S_densCOM = np.zeros([len(vdosCOMnorm)])
S_densROT = np.zeros([len(vdosROTnorm)])

for i in range(1,len(vdosROTnorm)):
  S_densCOM[i] = COMEntropySolidQM(vdos[0,i],vdosCOMnorm[i],T) * Na / df
  S_densROT[i] = ROTEntropySolidQM(vdos[0,i],vdosROTnorm[i],T) * Na / df

#######################
# Save data to files: #
#######################

base_name = sys.argv[7]

outfile = open(base_name + "_Sdens.dat","w+")
for i in range(len(Scom)):
  outfile.write("%f\t%f\t%f\t%f\n" % (vdos[0,i], S_densCOM[i], S_densROT[i], S_densCOM[i]+S_densROT[i]))
outfile.close()

outfile = open(base_name + "_entropy.dat","w+")
for i in range(len(Scom)):
  outfile.write("%f\t%f\t%f\t%f\n" % (vdos[0,i], Scom[i], Srot[i], Scom[i]+Srot[i]))
outfile.close()

outfile = open(base_name + "_vdos.dat","w+")
for i in range(len(Scom)):
  outfile.write("%f\t%f\t%f\t%f\n" % (vdos[0,i], vdosCOMnorm[i]*10**13, vdosROTnorm[i]*10**13, (vdosCOMnorm[i]+vdosROTnorm[i])*10**13))
outfile.close()

###################
# Plot entropies: #
###################

# error_tot = (1.0 - (Srot[-1] + Scom[-1])/ref_tot) * 100.0

# fig, ax = plt.subplots(figsize=(4, 3))
# x1, y1 = [0, 1200], [ref_tot, ref_tot]
# x2, y2 = [0, 1200], [ref_rot , ref_rot]
# x3, y3 = [0, 1200], [ref_trans, ref_trans]
# ax.plot(vdos[0], Srot,color='r',label='SRE rot')
# ax.plot(vdos[0], Scom,color='b',label='SRE trans')
# ax.plot(vdos[0], Scom+Srot,color='gray',label='SRE tot')
# ax.plot(x1,y1,color='gray',linestyle='--',label='total bulk (2PT)')
# ax.plot(x2,y2,color='r',linestyle='--',label='rot bulk (2PT)')
# ax.plot(x3,y3,color='b',linestyle='--',label='trans bulk (2PT)')
# ax.set(xlabel=r'Frequency ($\mathrm{cm^{-1}}$)', ylabel=r'Entropy ($\mathrm{J/K\cdot mol}$)')
# ax.text(100,12,'deviation from bulk: '+ str(int(ref_tot - (Srot[-1] + Scom[-1]))) + r'$\mathrm{J/K\cdot mol}$')
# plt.legend(loc='center right',fontsize=7)
# plt.xlim(0, 1200)
# plt.ylim(0, 60)
# plt.tight_layout()
# plt.savefig(base_name + "_entropy_SRE" + ".pdf")
# plt.show()

#########################
# Plot normalized VDOS: #
#########################

# fig, ax = plt.subplots(figsize=(4, 3))
# ax.plot(vdos[0],(vdosROTnorm+vdosCOMnorm)*10**13, color='black', label='total')
# ax.plot(vdos[0], vdosROTnorm*10**13, color='red', label='rotation')
# ax.plot(vdos[0], vdosCOMnorm*10**13, color='blue', label='translation')
# ax.set(xlabel=r'Frequency ($\mathrm{cm^{-1}}$)', ylabel=r'Spectral density (arb. u.)')
# plt.legend(loc='upper right',fontsize=9)
# plt.xlim(0, 1200)
# plt.ylim(0, 10)
# plt.tight_layout()
# plt.savefig(base_name + "_vdos" + ".pdf")
# plt.show()