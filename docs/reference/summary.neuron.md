# Summarize a neuron object

Builds a fuller summary of a neuron object (class `"neuron"`, as created
by
[`new.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/new.neuron.md)/[`load.rasters.as.neurons`](https://michaelbarkasi.github.io/neuronsDG/reference/load.rasters.as.neurons.md))
than
[`print.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/print.neuron.md)
gives: identifying metadata, whether trial data is loaded (and, if so,
its mean spike rate), whether autocorrelation has been computed, and
whether an exponential-decay fit (and its parameters) is available.

## Usage

``` r
# S3 method for class 'neuron'
summary(object, ...)
```

## Arguments

- object:

  A neuron object.

- ...:

  Unused; present for S3 consistency.

## Value

A `summary.neuron` object, printed by its own
[`print.summary.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/print.summary.neuron.md)
method.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
summary(nrn)
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
