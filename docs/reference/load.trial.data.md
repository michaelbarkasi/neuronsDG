# Load a trial data matrix into a neuron object

Loads a trial data matrix directly into a neuron object, bypassing
raster parsing. Mainly useful for constructing a neuron that shares
another neuron's already-loaded trial data.

## Usage

``` r
load.trial.data(nrn, trial_data)
```

## Arguments

- nrn:

  A neuron object.

- trial_data:

  A numeric matrix of trial data, rows as time bins and columns as
  trials (e.g. from
  [`fetch.trial.data`](https://michaelbarkasi.github.io/neuronsDG/reference/fetch.trial.data.md)).

## Value

`nrn`, invisibly. Modifies `nrn` in place.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
nrn_copy <- new.neuron(t_per_bin = fetch.id.data(nrn)$t_per_bin)
load.trial.data(nrn_copy, fetch.trial.data(nrn))
```
