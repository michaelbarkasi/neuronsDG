
# Demo data download-and-cache helpers ################################################################################
#
# The real demo datasets used by the vignettes (kilosort4 output across 7
# recording sessions, and two large compact spike rasters) are too large to
# bundle in the package (together, ~89MB). Instead they are hosted as GitHub
# release assets and downloaded once into a per-user cache directory
# (tools::R_user_dir("neuronsDG", "cache")), from where they are reused on
# every subsequent call. If no internet connection is available, or the
# download otherwise fails, these functions return NULL (after a warning)
# rather than erroring, so callers (vignettes, examples, tests) can detect
# that and skip gracefully instead of failing a check.

.neuronsDG_demo_data_repo <- "michaelbarkasi/neuronsDG"
.neuronsDG_demo_data_tag  <- "demo-data-v1"

demo_cache_dir <- function() {
  d <- tools::R_user_dir("neuronsDG", "cache")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  d
}

demo_data_url <- function(filename) {
  paste0(
    "https://github.com/", .neuronsDG_demo_data_repo, "/releases/download/",
    .neuronsDG_demo_data_tag, "/", filename
  )
}

# Download `filename` into the cache dir if not already there; returns the
# cached path, or NULL on failure (leaving no partial file behind).
download_demo_asset <- function(filename) {
  dest <- file.path(demo_cache_dir(), filename)
  if (file.exists(dest)) return(dest)

  ok <- tryCatch({
    utils::download.file(
      demo_data_url(filename), destfile = dest, mode = "wb", quiet = TRUE
    )
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)

  if (!ok || !file.exists(dest)) {
    unlink(dest)
    return(NULL)
  }
  dest
}

#' Download the kilosort4 demo dataset
#'
#' The kilosort4 demo dataset (7 recording sessions) used by
#' \code{vignette("tutorial_tau_est_kilosort", package = "neuronsDG")} is
#' hosted as a GitHub release asset rather than bundled in the package, to
#' keep the installed package small. This downloads and unpacks it once into
#' a per-user cache directory (see \code{\link[tools]{R_user_dir}}) and
#' reuses that cache on every later call.
#'
#' @return Path to the cached \code{kilo4demo} directory (suitable for
#'   \code{preprocess.kilo4()}'s \code{recording.folder} argument), or
#'   \code{NULL} if the data could not be downloaded (e.g. no internet
#'   access), with a warning explaining why.
#' @export
get_kilo4demo_data <- function() {
  cache_dir <- demo_cache_dir()
  target_dir <- file.path(cache_dir, "kilo4demo")
  if (dir.exists(target_dir)) return(target_dir)

  zip_path <- download_demo_asset("kilo4demo.zip")
  if (is.null(zip_path)) {
    warning(
      "Could not download the kilosort4 demo dataset (no internet access, ",
      "or the GitHub release asset is unavailable). See ?get_kilo4demo_data.",
      call. = FALSE
    )
    return(NULL)
  }

  utils::unzip(zip_path, exdir = cache_dir)
  unlink(zip_path)
  if (!dir.exists(target_dir)) return(NULL)
  target_dir
}

#' Download the cross-correlation demo dataset
#'
#' The compact spike raster used by
#' \code{vignette("tutorial_synchrony_crosscorr", package = "neuronsDG")} is
#' hosted as a GitHub release asset rather than bundled in the package, to
#' keep the installed package small. This downloads it once into a per-user
#' cache directory (see \code{\link[tools]{R_user_dir}}) and reuses that
#' cache on every later call.
#'
#' @return A data frame (read directly from the cached, gzip-compressed
#'   file), or \code{NULL} if the data could not be downloaded (e.g. no
#'   internet access), with a warning explaining why.
#' @export
get_cross_corr_demo_data <- function() {
  gz_path <- download_demo_asset("cross_corr_demo_data.csv.gz")
  if (is.null(gz_path)) {
    warning(
      "Could not download the cross-correlation demo dataset (no internet ",
      "access, or the GitHub release asset is unavailable). ",
      "See ?get_cross_corr_demo_data.",
      call. = FALSE
    )
    return(NULL)
  }
  utils::read.csv(gzfile(gz_path))
}

#' Download the 2022 spike-raster demo dataset
#'
#' The compact spike raster used by
#' \code{vignette("tutorial_tau_est_DG", package = "neuronsDG")} is hosted as
#' a GitHub release asset rather than bundled in the package, to keep the
#' installed package small. This downloads it once into a per-user cache
#' directory (see \code{\link[tools]{R_user_dir}}) and reuses that cache on
#' every later call.
#'
#' @return A data frame (read directly from the cached, gzip-compressed
#'   file), or \code{NULL} if the data could not be downloaded (e.g. no
#'   internet access), with a warning explaining why.
#' @export
get_spike_rasters_2022_data <- function() {
  gz_path <- download_demo_asset("spike_rasters_2022data.csv.gz")
  if (is.null(gz_path)) {
    warning(
      "Could not download the 2022 spike-raster demo dataset (no internet ",
      "access, or the GitHub release asset is unavailable). ",
      "See ?get_spike_rasters_2022_data.",
      call. = FALSE
    )
    return(NULL)
  }
  utils::read.csv(gzfile(gz_path))
}
