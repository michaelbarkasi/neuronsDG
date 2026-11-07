# Print a neuron summary

Formats a `summary.neuron` object (as returned by
[`summary.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/summary.neuron.md))
for console display.

## Usage

``` r
# S3 method for class 'summary.neuron'
print(x, ...)
```

## Arguments

- x:

  A `summary.neuron` object, from
  [`summary.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/summary.neuron.md).

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
print(summary(nrn)) # print() dispatches to print.summary.neuron() automatically
#> Neuron #1 (not_provided)
#>   type: generic 
#>   simulated: FALSE 
#>   unit_time: ms 
#>   t_per_bin: 10 ms 
#>   sample_rate: 10000 Hz 
#>   trial data loaded: TRUE 
#>     lambda:     0.1 spike per ms 
#>     lambda_bin: 1 spike per bin
#>   autocorrelation computed: FALSE 
#>   exponential decay fit: FALSE 
```
