# Add new cell type

This function adds a user-defined cell type to the current session. It's
just a wrapper for the Rcpp-exported `add_cell_type` function.
Technically, `cell_type` is a `struc` defined in the Rcpp backend of the
neurons package. They are essentially labeled lists with the following
entries: `type_name`, `valence`, `temporal_modulation_bias`,
`temporal_modulation_timeconstant`, `temporal_modulation_amplitude`,
`transmission_velocity`, `v_bound`, `dHdv_bound`, `I_spike`,
`coupling_scaling_factor`, `spike_potential`, `resting_potential`, and
`threshold`. Each session stores cell types in the Rcpp backend in an
`unordered_map` with `string` labels. All parameters come with
biologically realistic (and mathematically workable) default values,
except for `type_name` and `valence`.

## Usage

``` r
add.cell.type(
  type_name,
  valence,
  temporal_modulation_bias = 0.001,
  temporal_modulation_timeconstant = 1,
  temporal_modulation_amplitude = 0.005,
  transmission_velocity = 30000,
  v_bound = 85,
  dHdv_bound = 1.05e-06,
  I_spike = 1e-06,
  coupling_scaling_factor = 1e-07,
  spike_potential = 35,
  resting_potential = -70,
  threshold = -55
)
```

## Arguments

- type_name:

  Character string giving name of the cell type, e.g. "excitatory",
  "inhibitory", "PV", "SST", etc.

- valence:

  Valence of each neuron type, +1 for excitatory, -1 for inhibitory

- temporal_modulation_bias:

  Temporal modulation time (in ms) bias for each neuron type. Default
  value is 1e-3.

- temporal_modulation_timeconstant:

  Temporal modulation time (in ms) step for each neuron type. Default
  value is 1e0.

- temporal_modulation_amplitude:

  Temporal modulation time (in ms) cutoff for each neuron type. Default
  value is 5e-3.

- transmission_velocity:

  Transmission velocity (in microns/ms) for each neuron type. Default
  value is 30e3.

- v_bound:

  Potential bound, such that -v_bound \<= v_traces \<= v_bound, in
  unit_potential (mV), for each neuron in the network, based on its
  type. Default value is 85.0.

- dHdv_bound:

  Bound on derivative of metabolic energy wrt potential, such that
  dHdv_bound \> abs(dHdv), in mA, for each neuron in the network, based
  on its type. Default value is 1.05e-6.

- I_spike:

  Spike current, in mA. Default value is 1e-6 (i.e., 1 nA).

- coupling_scaling_factor:

  Controls how energy used in synaptic transmission compares to that
  used in spiking. Default value is 1e-7, meaning that synaptic
  transmission uses 0.00001 percent of the energy used in spiking.

- spike_potential:

  Magnitude of each spike, in mV. Default value is 35.0.

- resting_potential:

  Resting potential, in mV. Default value is -70.0.

- threshold:

  Spike threshold, in mV. Default value is -55.0.

## Value

Nothing.
