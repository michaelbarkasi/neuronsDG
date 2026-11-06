# Fetch a neuron's spike raster

Fetch a neuron's spike raster

## Usage

``` r
fetch.spike.raster(nrn)
```

## Arguments

- nrn:

  A neuron object.

## Value

A numeric matrix giving the neuron's spike raster (assumes it has
already been loaded).

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
fetch.spike.raster(nrn)
#>      time trial
#> [1,]   10     1
#> [2,]   20     1
#> [3,]   30     1
#> [4,]   40     1
#> [5,]   15     2
#> [6,]   25     2
#> [7,]   35     2
#> [8,]   45     2
```
