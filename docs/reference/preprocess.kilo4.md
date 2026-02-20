# Preprocess kilosort4 data

Function to preprocess kilosort4 data for use with neurons package
functions.

## Usage

``` r
preprocess.kilo4(
  trial_time_start = -100,
  trial_time_end = 2020,
  recording.folder = "data",
  meta_data = NULL,
  max_spikes = Inf,
  min_spikes = 0,
  min_trials = 0,
  pure_trials_only = TRUE,
  good_cells_only = TRUE,
  stim_responsive_only = TRUE,
  verbose = TRUE
)
```

## Arguments

- trial_time_start:

  Time to begin trial (ms), relative to stimulus onset

- trial_time_end:

  Time to end trial (ms), relative to stimulus onset

- recording.folder:

  List of paths to the kilosort4 output folders

- meta_data:

  Data frame with metadata for each recording (e.g., genotype,
  hemisphere), one recording per row; row names should match recording
  names and all columns should be covariates for later analysis

- max_spikes:

  Maximum number of total spikes for a cell to be kept

- min_spikes:

  Minimum number of total spikes for a cell to be kept

- min_trials:

  Minimum number of trials for a cell to be kept

- pure_trials_only:

  Keep only trials with no overlap?

- good_cells_only:

  Keep only cells marked as "good" in the cluster_group.tsv file?

- stim_responsive_only:

  Keep only cells marked as stimulus-responsive in the cluster_group.tsv
  file?

- verbose:

  Print report on imported data?

## Value

A list with three elements:

- spikes:

  data frame with one row per spike

- timeXtrial:

  list of matrices with rows as sample times and columns as trials, each
  element a zero if no spike at that time in that trial, and a one if a
  spike

- cluster.key:

  data frame with one row per neuron
