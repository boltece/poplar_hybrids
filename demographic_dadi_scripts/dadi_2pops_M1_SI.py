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


os.chdir('./easySFS/output_coastT_balsam/dadi/')
data = dadi.Spectrum.from_file("balsam-coastal_tricho.sfs")

ns=data.sample_sizes


## Strict Isolation model (M1)

pts_1 = [250,260,270]

def SI(params, ns, pts):
	nuB, nuT, T1 = params
	"""
	2-populations, strict isolation model with no gene flow.

	nuB:  Current size of balsam population, after split.
	nuT:  Current size of coastal_tricho population.
	T1:   Divergence time for split between coastal_tricho and balsam species
	ns:   Size of fs to generate.
	pts:  Number of points to use in grid for evaluation.
	"""
	xx = dadi.Numerics.default_grid(pts)
	phi = dadi.PhiManip.phi_1D(xx)
	phi = dadi.PhiManip.phi_1D_to_2D(xx,phi)
	phi = dadi.Integration.two_pops(phi,xx,T1,nu1=nuB, nu2=nuT)
	fs = dadi.Spectrum.from_phi(phi,ns,(xx,xx))
	return fs

func = SI



upper_bound = [10, 10, 10]
lower_bound = [1e-4,1e-4, 1e-4]
p0 = array([0.1,0.1, 0.1])


p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
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

matplotlib.pyplot.savefig("SI_unfolded_run1.pdf")

