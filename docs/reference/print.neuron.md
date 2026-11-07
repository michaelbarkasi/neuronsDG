# Print a neuron object

Terse console summary of a neuron object (class `"neuron"`, as created
by
[`new.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/new.neuron.md)/[`load.rasters.as.neurons`](https://michaelbarkasi.github.io/neuronsDG/reference/load.rasters.as.neurons.md)):
its identifying metadata, whether trial data has been loaded, and
whether autocorrelation has been computed. A neuron object is a plain S3
wrapper around the package's internal C++ representation, so
auto-printing it at the console (typing it by itself, rather than
calling [`print()`](https://rdrr.io/r/base/print.html) explicitly)
already goes through ordinary S3 dispatch to this method, with no
special handling needed.

## Usage

``` r
# S3 method for class 'neuron'
print(x, ...)
```

## Arguments

- x:

  A neuron object.

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
print(nrn) # or just: nrn
#> <neuron> #1 (not_provided) 
#>   type: generic 
#>   trial data loaded: TRUE 
#>     lambda:     0.1 spike per ms 
#>     lambda_bin: 1 spike per bin
#>   autocorrelation computed: FALSE 
```
