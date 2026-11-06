# Download the kilosort4 demo dataset

The kilosort4 demo dataset (7 recording sessions) used by
`vignette("tutorial_tau_est_kilosort", package = "neuronsDG")` is hosted
as a GitHub release asset rather than bundled in the package, to keep
the installed package small. This downloads and unpacks it once into a
per-user cache directory (see
[`R_user_dir`](https://rdrr.io/r/tools/userdir.html)) and reuses that
cache on every later call.

## Usage

``` r
get_kilo4demo_data()
```

## Value

Path to the cached `kilo4demo` directory (suitable for
[`preprocess.kilo4()`](https://michaelbarkasi.github.io/neuronsDG/reference/preprocess.kilo4.md)'s
`recording.folder` argument), or `NULL` if the data could not be
downloaded (e.g. no internet access), with a warning explaining why.

## Examples

``` r
if (FALSE) { # \dontrun{
kilo4_path <- get_kilo4demo_data()
list.files(kilo4_path)
} # }
```
