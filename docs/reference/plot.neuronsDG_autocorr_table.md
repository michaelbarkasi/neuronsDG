# Plot an autocorrelation summary table

`plot` method for the data frame returned by
[`process.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/process.autocorr.md).
Checks the assumption behind the sigma-fitting equation used to convert
spiking autocorrelation to a dichotomized-Gaussian sigma: the approach
is valid if autocorrelation (red lines) stays below mean firing rate
(blue line).

## Usage

``` r
# S3 method for class 'neuronsDG_autocorr_table'
plot(x, ...)
```

## Arguments

- x:

  A `neuronsDG_autocorr_table` object, i.e. the output of
  [`process.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/process.autocorr.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly. Produces a base-R plot as a side effect.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1,1,1,1, 1,1,1,1, 2,2,2,2),
  time_in_ms = c(10,20,30,40, 15,25,35,45, 12,22,32,42),
  trial      = c(1,1,1,1, 2,2,2,2, 1,1,1,1)
)
neurons <- load.rasters.as.neurons(raster_df, bin_size = 10)
batch <- process.autocorr(neurons)
plot(batch)
```
