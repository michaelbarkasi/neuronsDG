# Function to plot autocorrelation of neuron

Wrapper function to fetch and plot autocorrelation information from a
neuron object. Assumes the autocorrelation has already been computed.

## Usage

``` r
plot.autocorrelation(
 nrn, 
 plot_title = "Est. autocorr", 
 bias_term = 0, 
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
