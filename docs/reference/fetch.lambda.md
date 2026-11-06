# Fetch a neuron's mean spike rate

Fetch a neuron's mean spike rate

## Usage

``` r
fetch.lambda(nrn)
```

## Arguments

- nrn:

  A neuron object (assumes trial data has already been loaded).

## Value

A named numeric vector with elements `lambda` (mean spike rate per raw
time unit) and `lambda_bin` (mean spike rate per bin).

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
fetch.lambda(nrn)
#>     lambda lambda_bin 
#>        0.1        1.0 
```
