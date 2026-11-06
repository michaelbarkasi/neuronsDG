# Plot autocorrelation of a neuron

Fetches and plots autocorrelation information from a neuron object.
Assumes the autocorrelation has already been computed. Called by
`plot.neuron(nrn, type = "autocorrelation")`.

## Usage

``` r
plot_autocorrelation(
  nrn,
  plot_title = "Est. autocorr",
  bias_term = NULL,
  plot_time_cutoff = Inf,
  return_plot = FALSE
)
```

## Arguments

- nrn:

  Neuron object for which to plot autocorrelation.

- plot_title:

  Title for the plot (default: "Est. autocorr").

- bias_term:

  Bias term to plot as a horizontal line (default: NULL).

- plot_time_cutoff:

  Maximum lag (in bins) to display on the x-axis (default: Inf).

- return_plot:

  Logical indicating whether to return the ggplot object (TRUE) or print
  it (FALSE) (default: FALSE).

## Value

A ggplot object showing the estimated and fitted autocorrelation.

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
plot_autocorrelation(nrn)
```
