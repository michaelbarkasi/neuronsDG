# Download the 2022 spike-raster demo dataset

The compact spike raster used by
`vignette("tutorial_tau_est_DG", package = "neuronsDG")` is hosted as a
GitHub release asset rather than bundled in the package, to keep the
installed package small. This downloads it once into a per-user cache
directory (see [`R_user_dir`](https://rdrr.io/r/tools/userdir.html)) and
reuses that cache on every later call.

## Usage

``` r
get_spike_rasters_2022_data()
```

## Value

A data frame (read directly from the cached, gzip-compressed file), or
`NULL` if the data could not be downloaded (e.g. no internet access),
with a warning explaining why.

## Examples

``` r
if (FALSE) { # \dontrun{
spike_rasters <- get_spike_rasters_2022_data()
head(spike_rasters)
} # }
```
