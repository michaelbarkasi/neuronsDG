# Helper function for testing assumption of autocorrelation processing

The approach we're using is valid if autocorr (red lines) is below mean
firing rate (blue line). If so, then the equation we're using to convert
spiking autocorr to guassian sigma should work.

## Usage

``` r
test.sigma.assumption(autocor_results)
```

## Arguments

- autocor_results:

  Data frame of autocorrelation results from the `process.autocorr`
  function.

## Value

Nothing; this function produces a base-R plot for visual inspection.
