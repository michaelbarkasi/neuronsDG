# Modify existing cell type

This function modifies parameters of an existing cell type in the
current session. Parameters can be updated selectively. If the parameter
is not specified at all or is specified as `NULL`, the existing
parameter will be left in place.

## Usage

``` r
modify.cell.type(
  type_name,
  valence = NULL,
  temporal_modulation_bias = NULL,
  temporal_modulation_timeconstant = NULL,
  temporal_modulation_amplitude = NULL,
  transmission_velocity = NULL,
  v_bound = NULL,
  dHdv_bound = NULL,
  I_spike = NULL,
  coupling_scaling_factor = NULL,
  spike_potential = NULL,
  resting_potential = NULL,
  threshold = NULL
)
```

## Arguments

- type_name:

  Character string giving name of the cell type, e.g. "excitatory",
  "inhibitory", "PV", "SST", etc.

- valence:

  Valence of each neuron type, +1 for excitatory, -1 for inhibitory

- temporal_modulation_bias:

  Temporal modulation time (in ms) bias for each neuron type

- temporal_modulation_timeconstant:

  Temporal modulation time (in ms) step for each neuron type

- temporal_modulation_amplitude:

  Temporal modulation time (in ms) cutoff for each neuron type

- transmission_velocity:

  Transmission velocity (in microns/ms) for each neuron type

- v_bound:

  Potential bound, such that -v_bound \<= v_traces \<= v_bound, in mV,
  for each neuron in the network, based on its type

- dHdv_bound:

  Bound on derivative of metabolic energy wrt potential, such that
  dHdv_bound \> abs(dHdv), in mA, for each neuron in the network, based
  on its type

- I_spike:

  Spike current, in mA

- coupling_scaling_factor:

  Controls how energy used in synaptic transmission compares to that
  used in spiking

- spike_potential:

  Magnitude of each spike, in mV

- resting_potential:

  Resting potential, in mV

- threshold:

  Spike threshold, in mV

## Value

Nothing.
