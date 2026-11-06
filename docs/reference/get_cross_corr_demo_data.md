# Download the cross-correlation demo dataset

The compact spike raster used by
`vignette("tutorial_synchrony_crosscorr", package = "neuronsDG")` is
hosted as a GitHub release asset rather than bundled in the package, to
keep the installed package small. This downloads it once into a per-user
cache directory (see
[`R_user_dir`](https://rdrr.io/r/tools/userdir.html)) and reuses that
cache on every later call.

## Usage

``` r
get_cross_corr_demo_data()
```

## Value

A data frame (read directly from the cached, gzip-compressed file), or
`NULL` if the data could not be downloaded (e.g. no internet access),
with a warning explaining why.

## Examples

``` r
if (FALSE) { # \dontrun{
crosscorr_raster <- get_cross_corr_demo_data()
head(crosscorr_raster)
} # }
```
