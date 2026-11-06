# Plot spike raster of a neuron

Fetches and plots the spike raster from a neuron object. Assumes the
spike raster has already been loaded. Called by
`plot.neuron(nrn, type = "raster")` (the default plot for a neuron).

## Usage

``` r
plot_raster(nrn, plot_title = "Spike raster", zero_as_onset = TRUE, return_plot = FALSE)
```

## Arguments

- nrn:

  Neuron object for which to plot spike raster.

- plot_title:

  Title for the plot (default: "Spike raster").

- zero_as_onset:

  Logical indicating whether to plot a vertical line at time zero to
  indicate stimulus onset (default: TRUE).

- return_plot:

  Logical indicating whether to return the ggplot object (TRUE) or print
  it (FALSE) (default: FALSE).

## Value

A ggplot object showing the spike raster.

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1, 1, 1, 1, 1, 1, 1, 1),
  time_in_ms = c(10, 20, 30, 40, 15, 25, 35, 45),
  trial      = c(1, 1, 1, 1, 2, 2, 2, 2)
)
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]
plot_raster(nrn)
```
