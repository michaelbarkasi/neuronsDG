# Compute autocorrelation of single neuron object

This function computes the autocorrelation of a single neuron object
using specified parameters. It's a wrapper around the
`compute_autocorrelation` method of the neuron class.

## Usage

``` r
compute.autocorr(nrn, bin_count_action = "sum", max_lag = 0, use_raw = TRUE)
```

## Arguments

- nrn:

  Neuron object for which to compute autocorrelation.

- bin_count_action:

  Method for counting spikes in each bin when computing autocorrelation;
  one of "boolean", "mean", or "sum" (default: "sum").

- max_lag:

  Maximum lag (in units of the trial data) to compute autocorrelation;
  if 0, uses the number of time bins in the neuron's trial data
  (default: 0).

- use_raw:

  Logical indicating whether to use raw autocorrelation (true) or
  standard centered and normalized correlation (false) (default: TRUE).

## Value

None. Modifies neuron in place.
