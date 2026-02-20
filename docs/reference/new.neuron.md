# Initialize neuron

This function initializes a new neuron object with specified parameters.
The neuron objects are the main tools for running correlation analyses
of spike data.

## Usage

``` r
new.neuron(
  id_num = 0,
  recording_name = "not_provided",
  type = "generic",
  genotype = "not_provided",
  sex = "not_provided",
  hemi = "not_provided",
  region = "not_provided",
  age = "not_provided",
  sim = FALSE,
  unit_time = "ms",
  unit_sample_rate = "Hz",
  unit_data = "mV",
  t_per_bin = 10,
  sample_rate = 10000
)
```

## Arguments

- id_num:

  Numeric identifier for the neuron (default: 0).

- recording_name:

  Recording (if any) on which this neuron is based (default:
  "not_provided").

- type:

  Modeled electrophysiology of neuron, e.g. "generic", "blackbox" "LIF",
  "McCullochPitts", "excitatory", "inhibitory", etc. (default:
  "generic").

- genotype:

  Genotype of animal, e.g. "WT", "KO", "MECP2", "transgenic", etc.
  (default: "not_provided").

- sex:

  Sex of animal (default: "not_provided").

- hemi:

  Hemisphere of neuron in which neuron is located, e.g. "left", "right"
  (default: "not_provided").

- region:

  Brain region in which neuron is located, e.g. "V1", "M1", "CA1",
  "PFC", etc. (default: "not_provided").

- age:

  Age of animal, e.g. "P0", "P7", "P14", "adult", etc. (default:
  "not_provided").

- sim:

  Logical indicating if neuron is simulated (TRUE) or recorded (FALSE)
  (default: FALSE).

- unit_time:

  Unit of time for spike raster or other recording data, e.g. "ms", "s",
  etc. (default: "ms").

- unit_sample_rate:

  Unit of sample rate for spike raster or other recording data, e.g.
  "Hz", "kHz", etc. (default: "Hz").

- unit_data:

  Unit of data for spike raster or other recording data, e.g. "mV",
  "spike", etc. (default: "mV").

- t_per_bin:

  Time (in above units) per bin, e.g., 1 ms per bin (default: 10.0).

- sample_rate:

  Sample rate (in above units), e.g., 1e4 Hz (default: 1e4).

## Value

A new neuron object.
