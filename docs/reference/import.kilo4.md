# Import raw kilosort4 data

Base function to import raw kilosort4 and accompanying stimulus data.
The function looks for a csv file stim_data_file and looks for the
following files in folder_path/kilosort4:

- spike_positions.npy:

  2D array giving the x and y position of each spike

- spike_clusters.npy:

  integer giving the cluster number of each spike

- spike_times.npy:

  sample number the spike occurred at

- cluster_group.tsv:

  2D array giving status of each cluster (0=noise, 1=MUA, 2=Good,
  3=unsorted), hand-curated

- cluster_info.tsv:

  2D array giving the automatic output of kilosort4. This file not
  needed if cluster_group.tsv has data.

- includeVector.mat:

  MATLAB file giving whether each cluster is stimulus-responsive (1) or
  not (0)

## Usage

``` r
import.kilo4(
  folder_path,
  stim_data_file,
  recording_Fq_MATLAB = 30303,
  trial_time_start = -100,
  trial_time_end = 400,
  verbose = TRUE
)
```

## Arguments

- folder_path:

  Path to the kilosort4 output files

- stim_data_file:

  Path to the stimulus event data file

- recording_Fq_MATLAB:

  Sampling frequency of the MATLAB recording

- trial_time_start:

  Time to begin trial (ms), relative to stimulus onset

- trial_time_end:

  Time to end trial (ms), relative to stimulus onset

- verbose:

  Print report on imported data?

## Value

A list with two data frames: spikes (rows are spikes) and stim (rows are
stimulus onset times)
