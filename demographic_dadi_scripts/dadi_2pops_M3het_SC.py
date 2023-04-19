#!/usr/bin/env python
# coding: utf-8


import os
import matplotlib
matplotlib.use('PDF')
import dadi
import pylab
import matplotlib.pyplot as plt
import numpy as np
from numpy import array
from dadi import Misc,Spectrum,Numerics,PhiManip,Integration



os.chdir('/storage/work/ceb6313/easySFS/output_coastT_balsam/dadi/')
data = dadi.Spectrum.from_file("balsam-coastal_tricho.sfs")


ns=data.sample_sizes


# # Secondary contact model (SC) with two categories of loci experiencing heterogeneous migration rates throughout the genome


pts_1 = [250,260,270]

def SChet(params, ns, pts):
	nuB, nuT, T1, T2, m12, m21, m12h, m21h, P = params
	"""
	2-populations, secondary contact model with no gene flow during T1.

	nuB:  Current size of balsam population, after split.
	nuT:  Current size of coastal tricho population.
	m12: gene flow rate from coastal tricho into balsam (if het, then genomic island param)
	m21: gene flow rate from balsam into coastal tricho (if het, then genomic island param)
	T1:   Time for divergence between balsam and tricho species with no gene flow
	T2: Time interval with gene flow
	P:  proportion of the genome evolving neutrally
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	### Calculate the neutral spectrum
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuB, nu2=nuT,m12=0, m21=0)
	phi = dadi.Integration.two_pops(phi,xx,T2,nu1=nuB, nu2=nuT,m12=m12, m21=m21)
	fsUn = dadi.Spectrum.from_phi(phi,ns,(xx,xx))

	### Calculate the genomic island spectrum
	phiI = dadi.PhiManip.phi_1D(xx)
	phiI = dadi.PhiManip.phi_1D_to_2D(xx,phiI)
	phiI = dadi.Integration.two_pops(phiI, xx,T1, nu1=nuB, nu2=nuT, m12=0,m21=0)
	phiI = dadi.Integration.two_pops(phiI, xx,T2, nu1=nuB, nu2=nuT, m12=m12h,m21=m21h)
	fsI = dadi.Spectrum.from_phi(phiI,ns,(xx,xx))
	
	### Sum the two spectra in proportions P and 1-P, *without* integrating the effect of misorientation of ancestor (we used a consensus assignment for ancestral state)
	fs = P*fsUn + (1-P)*fsI
	return fs

func = SChet



upper_bound = [10, 10, 0.1, 0.1, 100, 100, 100, 100, 0.9]
lower_bound = [1e-4,1e-4, 1e-6, 1e-6, 1e-3, 1e-3, 1e-3, 1e-3, 0.1]
p0 = array([0.1,0.1, 0.001, 0.001, 1, 1, 1, 1,0.4])




p0 = dadi.Misc.perturb_params(p0, upper_bound=upper_bound,lower_bound=lower_bound)
print('starting set: {0}'.format(p0))
func_ex = dadi.Numerics.make_extrap_log_func(func)
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_1,lower_bound=lower_bound,upper_bound=upper_bound,verbose=len(p0), maxiter=25)



print('Best-fit parameters: {0}'.format(popt))


model = func_ex(popt, ns, pts_1)
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))


theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))



plt=dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=2, pop_ids =("balsam","coastal_tricho"))

matplotlib.pyplot.savefig("SChet_unfolded.pdf")

