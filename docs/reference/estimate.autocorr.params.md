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
  use_raw = TRUE,
  seed = 12345
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

- max_lag:

  Maximum lag (in raw, unbinned units of each neuron's own trial data)
  to compute autocorrelation for both the original and each
  dichotomized-Gaussian-simulated spike raster; if 0, uses the shortest
  raw trial length across `neuron_list` (default: 0).

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

- seed:

  Integer seed for the dichotomized-Gaussian simulation draws (default:
  12345). These draws are made by a self-contained C++ random number
  generator, not R's own, so
  [`set.seed()`](https://rdrr.io/r/base/Random.html) has no effect on
  them – reproducibility is controlled entirely by this argument. Each
  neuron in `neuron_list` gets its own seed, offset from this one so
  that neurons don't share a draw sequence.

## Value

A `neuronsDG_estimates` object: a list containing a data frame of
autocorrelation parameter estimates (one row per simulation), a data
frame of neuron identifiers (one row per neuron), and the number of
simulations run per neuron. See
[`print.neuronsDG_estimates`](https://michaelbarkasi.github.io/neuronsDG/reference/print.neuronsDG_estimates.md).

## Examples

``` r
# A well-conditioned example needs more trials/spikes than process.autocorr()'s
# minimal examples, so the dichotomized-Gaussian fit underneath isn't degenerate.
set.seed(42)
make_cell <- function(cell_id, hemi, n_trials = 20, trial_dur = 100, n_spikes = 10) {
  times <- unlist(lapply(
    seq_len(n_trials), function(tr) sort(sample(seq_len(trial_dur), n_spikes))
  ))
  data.frame(
    cell = cell_id, time_in_ms = times,
    trial = rep(seq_len(n_trials), each = n_spikes), hemi = hemi
  )
}
raster_df <- rbind(make_cell(1, "left"), make_cell(2, "right"))
neurons <- load.rasters.as.neurons(raster_df, bin_size = 10)

ests <- estimate.autocorr.params(
  neuron_list = neurons, n_trials_per_sim = 10, n_sims_per_neurons = 2, seed = 1
)
ests
#> <neuronsDG_estimates>
#>   neurons:           2 
#>   sims per neuron:   2 
#>   total simulations: 4 
#>   mean tau:          1.16 
#>   mean A:            0.3406 
```
