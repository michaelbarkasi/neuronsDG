# Installation instructions

## Installing the package

``` r

install.packages("neuronsDG")
```

**Note:** neuronsDG is not yet an accepted CRAN package. Until it is,
the command above will not work; install from source instead, as
described below.

## Installing from source

For contributors, or until the package is accepted on CRAN:

1.  Clone and `cd` into this git repo
2.  Make `build_install.sh` executable by running
    `chmod +x build_install.sh` in bash terminal
3.  In terminal, run `./build_install.sh` to build and install the
    package
4.  To ensure a clean start, run:

``` bash
rm -f src/*.o src/*.so
rm -rf neuronsDG.Rcheck
rm -f neuronsDG_*.tar.gz
./build_install.sh
```

## Build documentation

In R, run:

``` r

devtools::document()
pkgdown::build_site()
pkgdown::build_article("<article name>")
```

**Note:** these functions call out to [pandoc](https://pandoc.org/) to
render each article’s final HTML. RStudio bundles its own copy of pandoc
and finds it automatically, so the above just works in an RStudio
session. From a plain R or Rscript session outside RStudio, pandoc may
not be on the system `PATH`; if
[`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
fails with “Pandoc not available”, point R at a pandoc installation
first, e.g. (macOS, using the copy bundled with RStudio):

``` r

Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64")
```

(use the `x86_64` copy instead of `aarch64` on an Intel Mac), or install
pandoc directly from [pandoc.org](https://pandoc.org/installing.html).
