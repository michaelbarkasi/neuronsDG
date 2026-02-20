# Load neuron data from raster file

This function loads spike raster data from a data frame or CSV file and
converts each unique cell into a neuron object. The raster data frame
must contain columns for cell identifier, spike time in milliseconds,
and trial number. Optional metadata columns can also be included.

## Usage

``` r
load.rasters.as.neurons(
  raster_df,
  bin_size = 10,
  min_duration = 0,
  max_displacement = 1e+09,
  sample_rt = 10000,
  time_cutoff = Inf
)
```

## Arguments

- raster_df:

  Data frame (or file name to csv importable as such), each row a spike;
  must have columns: cell, time_in_ms, trial; optional columns:
  recording_name, hemisphere, genotype, sex, region, age.

- bin_size:

  Size of time bins in milliseconds (default: 10).

- min_duration:

  Minimum trial duration in milliseconds (default: 0).

- max_displacement:

  Max displacement of trial center in milliseconds (default: 1e9, i.e.,
  "infinity").

- sample_rt:

  Sample rate in the default unit for neuron objects, Hz (default: 1e4).

- time_cutoff:

  Maximum time (in ms) to include spikes; spikes occurring after this
  time will be excluded (default: Inf).

## Value

A list of neuron objects, one per unique cell in the raster data.
