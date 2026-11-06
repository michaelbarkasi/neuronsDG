# S3/S4 methods for neuron objects

`print`, `summary`, and `plot` methods for neuron objects (class
`"Rcpp_neuron"`, as created by
[`new.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/new.neuron.md)/[`load.rasters.as.neurons`](https://michaelbarkasi.github.io/neuronsDG/reference/load.rasters.as.neurons.md)).
`plot.Rcpp_neuron` dispatches to
[`plot_raster`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_raster.md)
or
[`plot_autocorrelation`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_autocorrelation.md)
depending on `type`, defaulting to the raster view. A neuron object is
formally S4, so auto-printing it at the console (typing it by itself,
rather than calling [`print()`](https://rdrr.io/r/base/print.html)
explicitly) goes through an S4 `show` method instead, which simply calls
`print.Rcpp_neuron`.

## Usage

``` r
# S3 method for class 'Rcpp_neuron'
print(x, ...)

# S3 method for class 'Rcpp_neuron'
summary(object, ...)

# S3 method for class 'summary.neuron'
print(x, ...)

# S3 method for class 'Rcpp_neuron'
plot(x, type = c("raster", "autocorrelation"), ...)
```

## Arguments

- x:

  A neuron object (or, for `print.summary.neuron`, a `summary.neuron`
  object from `summary.Rcpp_neuron`).

- ...:

  Unused for `print`/`summary`; present for S3 consistency. For `plot`,
  additional arguments passed on to
  [`plot_raster`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_raster.md)
  or
  [`plot_autocorrelation`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_autocorrelation.md).

- object:

  A neuron object.

- type:

  One of "raster" (default) or "autocorrelation", specifying which view
  of the neuron to plot.

## Value

`print.Rcpp_neuron` and `print.summary.neuron` return their argument,
invisibly. `summary.Rcpp_neuron` returns a `summary.neuron` object,
printed by its own `print` method. `plot.Rcpp_neuron` returns whatever
[`plot_raster`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_raster.md)/[`plot_autocorrelation`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_autocorrelation.md)
returns (a ggplot object if `return_plot = TRUE`, otherwise `NULL`,
invisibly).

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

print(nrn)
#> <neuron> #1 (not_provided) 
#>   type: generic 
#>   trial data loaded: TRUE 
#>     lambda:     0.1 spike per ms 
#>     lambda_bin: 1 spike per bin
#>   autocorrelation computed: TRUE 
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
#>   autocorrelation computed: TRUE 
#>   exponential decay fit: TRUE 
#>     A:         0.9241 
#>     tau:       1.565 ms 
#>     bias_term: 1 
plot(nrn) # raster view (the default)

plot(nrn, type = "autocorrelation")
```
