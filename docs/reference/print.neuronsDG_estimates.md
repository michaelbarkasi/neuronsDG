# Print an autocorrelation-parameter estimate set

Print method for the list returned by
[`estimate.autocorr.params`](https://michaelbarkasi.github.io/neuronsDG/reference/estimate.autocorr.params.md).

## Usage

``` r
# S3 method for class 'neuronsDG_estimates'
print(x, ...)
```

## Arguments

- x:

  A `neuronsDG_estimates` object, from
  [`estimate.autocorr.params`](https://michaelbarkasi.github.io/neuronsDG/reference/estimate.autocorr.params.md).

- ...:

  Unused; present for S3 consistency.

## Value

`x`, invisibly.

## Examples

``` r
set.seed(42)
times <- unlist(lapply(1:20, function(tr) sort(sample(1:100, 10))))
raster_df <- data.frame(cell = 1, time_in_ms = times, trial = rep(1:20, each = 10))
nrn <- load.rasters.as.neurons(raster_df, bin_size = 10)[["neuron_1"]]

ests <- estimate.autocorr.params(
  neuron_list = list(nrn), n_trials_per_sim = 10, n_sims_per_neurons = 2, seed = 1
)
print(ests)
#> <neuronsDG_estimates>
#>   neurons:           1 
#>   sims per neuron:   2 
#>   total simulations: 2 
#>   mean tau:          1.16 
#>   mean A:            0.3406 
```
