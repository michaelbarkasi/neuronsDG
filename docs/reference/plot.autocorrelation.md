# Function to plot autocorrelation of neuron

Wrapper function to fetch and plot autocorrelation information from a
neuron object. Assumes the autocorrelation has already been computed.

## Usage

``` r
# S3 method for class 'autocorrelation'
plot(nrn, plot_title = "Est. autocorr", bias_term = 0, plot_time_cutoff = Inf)
```

## Arguments

- nrn:

  Neuron object for which to plot autocorrelation.

- plot_title:

  Title for the plot (default: "Est. autocorr").

- bias_term:

  Bias term to plot as a horizontal line (default: 0).

- plot_time_cutoff:

  Maximum lag (in bins) to display on the x-axis (default: Inf).

## Value

A ggplot object showing the estimated and fitted autocorrelation.
