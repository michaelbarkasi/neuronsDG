# Function to analyze autocorrelation estimates with bootstrapping

This function takes the output from `estimate.autocorr.params` and
performs a bootstrap analysis to estimate the distribution of the mean
time constant (tau) for different levels of specified covariates.

## Usage

``` r
analyze.autocorr(ests, covariate, n_bs = 10000)
```

## Arguments

- ests:

  A `neuronsDG_estimates` object (the output of
  `estimate.autocorr.params`), or a list of such objects.

- covariate:

  Character string or vector of character strings specifying the
  covariate(s) to analyze (e.g., "hemi").

- n_bs:

  Number of bootstrap resamples to perform (default: 1e4).

## Value

A `neuronsDG_analysis` object: a list with two elements:

- resamples:

  data frame with bootstrap resamples of the mean tau for each
  combination of covariate levels

- distribution_plot:

  ggplot object showing the distribution of mean tau for each
  combination of covariate levels

See
[`neuronsDG_analysis-methods`](https://michaelbarkasi.github.io/neuronsDG/reference/neuronsDG_analysis-methods.md)
for its `print`/`plot` methods.

## Examples

``` r
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
result <- analyze.autocorr(ests, covariate = "hemi", n_bs = 50)
result
#> <neuronsDG_analysis>
#>   bootstrap resamples: 50 
#>   covariate groups:    left, right 
#>     left                 mean tau = 1.16
#>     right                mean tau = 1.16
plot(result)
```
