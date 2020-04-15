# sre_entropy
This script computes the configurational entropy of water based on the spectral density.
The method used here is the Spectrally Resolved Estimation (SRE) of the entropy
as previously described by Schlitter & Massarczyk (http://arxiv.org/abs/1909.04726).
Usually, due to the diffusivity of bulk water a purely harmonic approach 
severely underestimates the entropy of water. However, in confined systems and
for strongly slowed-down hydration water, a harmonic approach provides a lower bound
that is sufficiently close to the true entropy.

Input files need to be spectral densities (DoS) computed from the (mass-weighted)
velocity auto-correlation function, i.e.
  - translational DoS (from mass-weighted cartesian velocities)
  - rotational DoS (individual 3 angular velocities along
    the molecular axes of water molecules.

Output includes the smoothed spectra, the entropy density per unit frequency as well as
the entropy build-up as a function of frequency.
