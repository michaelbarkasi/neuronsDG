# neuronsDG: Modelling neurons via dichotomized Gaussians

neuronsDG is an R/C++ package for simulating neurons, currently under
development at [Oviedo Lab](https://oviedolab.org/). The package is
built around the core object class **neuron**, a C++ class accessible
through R via Rcpp modules and wrappers. This class is meant for
modeling the spiking of single neurons with dichotomized Gaussians.

![Overlay of spike trains with autocorrelation and dichotomized
Gaussian](reference/figures/DG_autocorr.png)

Artistic overlay of spike trains with autocorrelation and a dichotomized
Gaussian with covariance. This package uses the latter to generate
simulations of the former.

Copyright (C) 2026, Michael Barkasi <barkasi@wustl.edu>
