
# Title Page ###########################################################################################################
# 
# Neurons: A framework for neuron modelling, simulation, and analysis
# 
# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2024
#   
# Front matter #########################################################################################################

# Prelim code:
#   (none)

# Prelim comments and warnings: 
#   1. This script defines the R6 object class Neuron and related helper R6 classes, along with methods for them. 
#         Neuron objects and their methods are intended to be a flexible, common framework for integrating, 
#         analyzing, and simulating different kinds of neuron data. 

#library(Rcpp)
#library(RcppEigen)

# Source Cpp ###########################################################################################################

#sourceCpp("src/neuron.cpp")

.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("neuron", TRUE)
}

new_neuron <- function(
  id_num = 0, 
  recording_name = "not_provided", 
  type = "generic", 
  hemi = "not_provided",
  sim = FALSE, 
  unit_time = "ms", 
  unit_sample_rate = "Hz", 
  unit_data = "mV", 
  t_per_bin = 1.0, 
  sample_rate = 1e4
) {
  neuron <- new(
    neuron, 
    id_num, recording_name, type, hemi, sim, unit_time, unit_sample_rate, unit_data, t_per_bin, sample_rate
  )
  return(neuron)
}
