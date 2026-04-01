# Function to estimate autocorrelation parameters using dichotomized Gaussian simulations

This function performs the same estimate-and-fit procedure as
`process.autocorr`, but does so multiple times for each neuron using
simulated spike rasters generated from a dichotomized Gaussian model of
that neuron.

## Usage

``` r
estimate.autocorr.params(
  neuron_list,
  n_trials_per_sim = 300,
  n_sims_per_neurons = 100,
  bin_count_action = "sum",
  max_lag = 0,
  A0 = 0.001,
  tau0 = 1,
  ctol = 1e-08,
  max_evals = 500,
  use_raw = TRUE
)
```

## Arguments

- neuron_list:

  An R list of neuron objects.

- n_trials_per_sim:

  Number of trials to simulate for each simulation (default: 300).

- n_sims_per_neurons:

  Number of simulations to run for each neuron (default: 100).

- bin_count_action:

  Method for counting spikes in each bin when computing autocorrelation;
  one of "boolean", "mean", or "sum" (default: "sum").

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

- use_raw:

  Logical indicating whether to use raw autocorrelation (true) or
  standard centered and normalized correlation (false) (default: TRUE).

## Value

A list containing a data frame of autocorrelation parameter estimates
(one row per simulation), a data frame of neuron identifiers (one row
per neuron), and the number of simulations run per neuron.
