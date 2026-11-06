# Function to process autocorrelation of neuron list

First estimates autocorrelation of each neuron in the list, then fits an
exponential decay function to the estimated autocorrelation, and finally
collects and returns summary statistics for each neuron.

## Usage

``` r
process.autocorr(
  neuron_list,
  bin_count_action = "sum",
  max_lag = 0,
  A0 = 0.001,
  tau0 = 1,
  ctol = 1e-08,
  max_evals = 500,
  check_autofiring_ratio = FALSE,
  print_plots = FALSE,
  plot_time_cutoff = Inf,
  use_raw = TRUE
)
```

## Arguments

- neuron_list:

  An R list of neuron objects.

- bin_count_action:

  Method for counting spikes in each bin when computing autocorrelation;
  one of "boolean", "mean", or "sum" (default: "sum").

- max_lag:

  Maximum lag (in units of the trial data) to compute autocorrelation;
  if 0, uses the number of time bins in the neuron's trial data
  (default: 0).

- A0:

  Initial guess for amplitude parameter of exponential decay function
  (default: 0.001).

- tau0:

  Initial guess for time constant parameter of exponential decay
  function (default: 1.0).

- ctol:

  Convergence tolerance for fitting exponential decay function (default:
  1e-8).

- max_evals:

  Maximum number of evaluations for fitting exponential decay function
  (default: 500).

- check_autofiring_ratio:

  Logical indicating whether to check the assumption that
  autocorrelation values are below the mean firing rate with
  `test.sigma.assumption` function (default: FALSE).

- print_plots:

  Logical indicating whether to print autocorrelation plots for each
  neuron (default: FALSE).

- plot_time_cutoff:

  Maximum lag (in bins) to display on the x-axis of plots (default:
  Inf).

- use_raw:

  Logical indicating whether to use raw autocorrelation (true) or
  standard centered and normalized correlation (false) (default: TRUE).

## Value

A `neuronsDG_autocorr_table` object: a data frame with one row per
neuron and columns for lambda (mean spike rate per ms and per bin),
amplitude (A), time constant (tau), bias term, first autocorrelation
value, maximum autocorrelation value, mean autocorrelation value, and
minimum autocorrelation value. See
[`plot.neuronsDG_autocorr_table`](https://michaelbarkasi.github.io/neuronsDG/reference/plot.neuronsDG_autocorr_table.md).

## Examples

``` r
raster_df <- data.frame(
  cell       = c(1,1,1,1, 1,1,1,1, 2,2,2,2),
  time_in_ms = c(10,20,30,40, 15,25,35,45, 12,22,32,42),
  trial      = c(1,1,1,1, 2,2,2,2, 1,1,1,1)
)
neurons <- load.rasters.as.neurons(raster_df, bin_size = 10)
process.autocorr(neurons)
#>       cell lambda_ms lambda_bin         A      tau bias_term autocorr1
#> 1 neuron_1 0.1000000   1.000000 0.4613628 1.682201  1.000000  1.461363
#> 2 neuron_2 0.1333333   1.333333 0.3286061 1.079589  1.777778  2.106384
#>   max_autocorr mean_autocorr min_autocorr
#> 1     1.001209      1.000606     1.000003
#> 2     1.777809      1.777793     1.777778
```
