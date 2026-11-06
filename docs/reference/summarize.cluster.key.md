# Summarize cluster key by covariates

Function to summarize the `cluster.key` output of `preprocess.kilo4` by
covariates.

## Usage

``` r
summarize.cluster.key(key, covariate_list)
```

## Arguments

- key:

  The `cluster.key` data frame output of `preprocess.kilo4`

- covariate_list:

  List of covariates to summarize

## Value

A data frame with one row per combination of covariates, and columns
giving the number of cells, mean number of spikes, and mean number of
trials for that combination

## Examples

``` r
if (FALSE) { # \dontrun{
kilo4_path <- get_kilo4demo_data()
preprocessed <- preprocess.kilo4(recording.folder = kilo4_path)
summarize.cluster.key(preprocessed$cluster.key, covariate_list = c("hemi"))
} # }
```
