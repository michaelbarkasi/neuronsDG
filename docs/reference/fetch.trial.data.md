# Fetch a neuron's trial data matrix

Fetch a neuron's trial data matrix

## Usage

``` r
fetch.trial.data(nrn)
```

## Arguments

- nrn:

  A neuron object.

## Value

A numeric matrix of trial data, rows as time bins and columns as trials.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
fetch.trial.data(nrn)
#>        trial1 trial2
#> time1       1      0
#> time2       0      0
#> time3       0      0
#> time4       0      0
#> time5       0      1
#> time6       0      0
#> time7       0      0
#> time8       0      0
#> time9       0      0
#> time10      1      0
#> time11      0      0
#> time12      0      0
#> time13      0      0
#> time14      0      0
#> time15      0      1
#> time16      0      0
#> time17      0      0
#> time18      0      0
#> time19      0      0
#> time20      1      0
#> time21      0      0
#> time22      0      0
#> time23      0      0
#> time24      0      0
#> time25      0      1
#> time26      0      0
#> time27      0      0
#> time28      0      0
#> time29      0      0
#> time30      1      0
#> time31      0      0
#> time32      0      0
#> time33      0      0
#> time34      0      0
#> time35      0      1
#> time36      0      0
#> time37      0      0
#> time38      0      0
#> time39      0      0
#> time40      0      0
```
