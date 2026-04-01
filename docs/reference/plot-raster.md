# Function to plot spike raster of neuron

Wrapper function to fetch and plot spike raster from a neuron object.
Assumes the spike raster has already been loaded.

## Usage

``` r
plot.raster(nrn, plot_title = "Spike raster", zero_as_onset = TRUE, return_plot = FALSE)
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
