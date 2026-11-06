# Helper function for testing assumption of autocorrelation processing

The approach we're using is valid if autocorr (red lines) is below mean
firing rate (blue line). If so, then the equation we're using to convert
spiking autocorr to guassian sigma should work. Thin wrapper around
[`plot.neuronsDG_autocorr_table`](https://michaelbarkasi.github.io/neuronsDG/reference/plot.neuronsDG_autocorr_table.md),
which does the actual plotting.

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

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1,1,1,1, 1,1,1,1, 2,2,2,2),
  time_in_ms = c(10,20,30,40, 15,25,35,45, 12,22,32,42),
  trial      = c(1,1,1,1, 2,2,2,2, 1,1,1,1)
)
neurons <- load.rasters.as.neurons(raster_df, bin_size = 10)
batch <- process.autocorr(neurons)
test.sigma.assumption(batch)
```
