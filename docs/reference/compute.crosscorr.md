# Compute cross-correlation between two neurons

Computes the cross-correlation (raw or Pearson) between a reference
neuron's spike train and a lagged copy of a comparison neuron's spike
train. It's a wrapper around the `compute_crosscorrelation_R` method of
the neuron class.

## Usage

``` r
compute.crosscorr(
  nrn,
  nrn_compare,
  bin_count_action = "sum",
  max_lag,
  use_raw = TRUE,
  verbose = FALSE
)
```

## Arguments

- nrn:

  Reference neuron object.

- nrn_compare:

  Comparison neuron object; must have trial length at least as long as
  `nrn`'s.

- bin_count_action:

  Method for counting spikes in each bin; one of "boolean", "mean", or
  "sum" (default: "sum").

- max_lag:

  Maximum lag (in units of the trial data) to compute cross-correlation
  for.

- use_raw:

  Logical indicating whether to use raw cross-correlation (true) or
  standard centered and normalized correlation (false) (default: TRUE).

- verbose:

  Logical indicating whether to print progress information (default:
  FALSE).

## Value

A numeric vector of cross-correlation values at lags from zero to
`max_lag`.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 2, 2, 2, 2),
  time_in_ms = c(10, 20, 30, 40, 12, 22, 32, 42),
  trial      = c(1, 1, 1, 1, 1, 1, 1, 1)
)
neurons <- load.rasters.as.neurons(raster_df, bin_size = 1)
compute.crosscorr(neurons[["neuron_1"]], neurons[["neuron_2"]], max_lag = 5)
#> [1] 0.09230769 0.00000000 0.00000000 0.00000000 0.00000000
```
