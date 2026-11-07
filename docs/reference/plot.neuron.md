# Plot a neuron object

Default `plot` method for neuron objects (class `"neuron"`, as created
by
[`new.neuron`](https://michaelbarkasi.github.io/neuronsDG/reference/new.neuron.md)/[`load.rasters.as.neurons`](https://michaelbarkasi.github.io/neuronsDG/reference/load.rasters.as.neurons.md)).
Dispatches to
[`plot_raster`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_raster.md)
or
[`plot_autocorrelation`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_autocorrelation.md)
depending on `type`, defaulting to the raster view.

## Usage

``` r
# S3 method for class 'neuron'
plot(x, type = c("raster", "autocorrelation"), ...)
```

## Arguments

- x:

  A neuron object.

- type:

  One of "raster" (default) or "autocorrelation", specifying which view
  of the neuron to plot.

- ...:

  Additional arguments passed on to
  [`plot_raster`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_raster.md)
  or
  [`plot_autocorrelation`](https://michaelbarkasi.github.io/neuronsDG/reference/plot_autocorrelation.md).

## Value

Whatever
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

plot(nrn) # raster view (the default)

plot(nrn, type = "autocorrelation")
```
