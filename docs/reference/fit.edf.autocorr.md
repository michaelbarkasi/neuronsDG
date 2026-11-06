# Fit exponential decay function to autocorrelation of single neuron object

This function fits an exponential decay function to the autocorrelation
of a single neuron object using specified parameters. It's a wrapper
around the `fit_autocorrelation` method of the neuron class.

## Usage

``` r
fit.edf.autocorr(nrn, A0 = 0.001, tau0 = 1, ctol = 1e-08, max_evals = 500)
```

## Arguments

- nrn:

  Neuron object for which to fit autocorrelation.

- A0:

  Initial guess for amplitude parameter of exponential decay function
  (default: 0.001).

- tau0:

  Initial guess for time constant parameter of exponential decay
  function (default: 1.0).

- ctol:

  Convergence tolerance for fitting exponential decay function (default:
  1e-8).

- max_evals:

  Maximum number of evaluations for fitting exponential decay function
  (default: 500).

## Value

None. Modifies neuron in place.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
compute.autocorr(nrn)
fit.edf.autocorr(nrn)
fetch.edf.parameters(nrn)
#>         A       tau bias_term 
#> 0.9241097 1.5651261 1.0000000 
```
