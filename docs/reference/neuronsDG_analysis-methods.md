# Methods for autocorrelation bootstrap analyses

`print` and `plot` methods for the list returned by
[`analyze.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/analyze.autocorr.md).
`plot.neuronsDG_analysis` simply displays the stored distribution plot.

## Usage

``` r
# S3 method for class 'neuronsDG_analysis'
print(x, ...)

# S3 method for class 'neuronsDG_analysis'
plot(x, ...)
```

## Arguments

- x:

  A `neuronsDG_analysis` object, from
  [`analyze.autocorr`](https://michaelbarkasi.github.io/neuronsDG/reference/analyze.autocorr.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly. `plot.neuronsDG_analysis` prints `x$distribution_plot`
as a side effect.

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
print(result)
#> <neuronsDG_analysis>
#>   bootstrap resamples: 50 
#>   covariate groups:    left, right 
#>     left                 mean tau = 1.16
#>     right                mean tau = 1.16
plot(result)
```
