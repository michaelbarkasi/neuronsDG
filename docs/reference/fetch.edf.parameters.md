# Fetch a neuron's fitted exponential-decay parameters

Fetch a neuron's fitted exponential-decay parameters

## Usage

``` r
fetch.edf.parameters(nrn)
```

## Arguments

- nrn:

  A neuron object.

## Value

A named numeric vector with elements `A`, `tau`, and `bias_term`
(assumes a fit has already been run with
[`fit.edf.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/fit.edf.autocorr.md)).

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
