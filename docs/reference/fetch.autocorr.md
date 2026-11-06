# Fetch a neuron's computed empirical autocorrelation

Fetch a neuron's computed empirical autocorrelation

## Usage

``` r
fetch.autocorr(nrn)
```

## Arguments

- nrn:

  A neuron object.

## Value

A numeric vector of autocorrelation values by lag (assumes
autocorrelation has already been computed with
[`compute.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/compute.autocorr.md)).

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
compute.autocorr(nrn)
fetch.autocorr(nrn)
#> [1] 1.2500000 1.0000000 0.8333333 0.6562500
```
