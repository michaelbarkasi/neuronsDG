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

  Output from `estimate.autocorr.params`, or a list of such outputs.

- covariate:

  Character string or vector of character strings specifying the
  covariate(s) to analyze (e.g., "hemi").

- n_bs:

  Number of bootstrap resamples to perform (default: 1e4).

## Value

A list with two elements:

- resamples:

  data frame with bootstrap resamples of the mean tau for each
  combination of covariate levels

- distribution_plot:

  ggplot object showing the distribution of mean tau for each
  combination of covariate levels
