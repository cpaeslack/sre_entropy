# sre_entropy
This script computes the configurational entropy of water based on the spectral density,
as previously described by Schlitter & Massarczyk (http://arxiv.org/abs/1909.04726).

Input files need to be spectral densities (DoS) computed from the velocity auto-correlation function, 
i.e.
  - translational DoS (from mass-weighted cartesian velocities)
  - rotational DoS (individual 3 angular velocities along
    the molecular axes of water molecules
