
# Introduction

Neurons is an R/C++ package for simulating neurons and networks of neurons, currently under development at [Oviedo Lab](https://oviedolab.org/). The package is built around a few core object classes, all of which are C++ class accessible through R via Rcpp modules and wrappers.

- **neuron**: As the name suggests, intended for modeling single neurons. 
- **circuit**: For simulating spiking neural networks with Growth-transform dynamical systems. (Under development)
- **motif**: For specifying general patterns of connectivity between neurons in a circuit, i.e., a "circuit motif" or "network topology". (Under development)

Copyright (C) 2025, Michael Barkasi
barkasi@wustl.edu