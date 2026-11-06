# Regression values below were captured by running these exact calls against
# inst/extdata/spike_rasters_2022data.csv (the data used by
# vignettes/tutorial_tau_est_DG.Rmd) with set.seed(12345), bin_size = 10.0, on
# the current package build. They pin CURRENT behavior, not a certified
# ground truth.

test_that("compute.autocorr() + fit.edf.autocorr() match pinned values for the high-autocorrelation demo cell", {
  # fit.edf.autocorr()'s NLopt fit on real, noisy spike data is not bit-reproducible
  # across compilers/math libraries (known, accepted behavior -- see neuron::fit_autocorrelation()
  # in src/neuron.cpp); skip on CRAN's diverse build machines rather than pin an exact value that
  # can legitimately differ by platform.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()
  cell_high <- neurons[[16]]

  compute.autocorr(cell_high, use_raw = TRUE)
  fit.edf.autocorr(cell_high)
  edf <- cell_high$fetch_EDF_parameters()

  expect_equal(unname(edf["A"]), 0.1269031385, tolerance = 1e-5)
  expect_equal(unname(edf["tau"]), 69.3777201221, tolerance = 1e-3)
  expect_equal(unname(edf["bias_term"]), 0.0110591807, tolerance = 1e-6)

  autocorr <- cell_high$fetch_autocorr_R()
  expect_length(autocorr, 154)
  expect_equal(
    head(autocorr, 5),
    c(0.1720678702, 0.1242008765, 0.1005253776, 0.0876020265, 0.0792730489),
    tolerance = 1e-6
  )
})

test_that("compute.autocorr() + fit.edf.autocorr() match pinned values for the low-autocorrelation demo cell", {
  # See the high-autocorrelation test above for why this is platform-sensitive and CRAN-skipped.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()
  cell_low <- neurons[[1]]

  compute.autocorr(cell_low, use_raw = TRUE)
  fit.edf.autocorr(cell_low)
  edf <- cell_low$fetch_EDF_parameters()

  expect_equal(unname(edf["A"]), 0.0123107279, tolerance = 1e-5)
  expect_equal(unname(edf["tau"]), 37.4334419605, tolerance = 1e-3)
  expect_equal(unname(edf["bias_term"]), 0.0004492891, tolerance = 1e-6)
})

test_that("the demo's 'high autocorrelation' cell has a larger tau than the 'low autocorrelation' cell", {
  # This is the qualitative claim vignettes/tutorial_tau_est_DG.Rmd makes about
  # cell 16 ('long horizontal streaks') vs cell 1 ('more random distribution')
  neurons <- dg_demo_neurons()
  compute.autocorr(neurons[[16]]); fit.edf.autocorr(neurons[[16]])
  compute.autocorr(neurons[[1]]);  fit.edf.autocorr(neurons[[1]])

  tau_high <- neurons[[16]]$fetch_EDF_parameters()["tau"]
  tau_low  <- neurons[[1]]$fetch_EDF_parameters()["tau"]
  expect_gt(unname(tau_high), unname(tau_low))
})

test_that("process.autocorr() matches pinned per-neuron values and expected shape (regression)", {
  # See the high-autocorrelation test above for why the pinned A/tau values below are
  # platform-sensitive and CRAN-skipped.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()
  result <- process.autocorr(neurons)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 41)
  expect_equal(
    colnames(result),
    c("cell", "lambda_ms", "lambda_bin", "A", "tau", "bias_term",
      "autocorr1", "max_autocorr", "mean_autocorr", "min_autocorr")
  )

  row16 <- result[result$cell == "neuron_16", ]
  expect_equal(row16$tau, 69.3777286065, tolerance = 1e-3)
  expect_equal(row16$A, 0.1269031299, tolerance = 1e-5)

  row1 <- result[result$cell == "neuron_1", ]
  expect_equal(row1$tau, 37.4363218591, tolerance = 1e-3)
  expect_equal(row1$A, 0.0123103612, tolerance = 1e-5)
})

test_that("process.autocorr() uses a shared max_lag across the whole neuron list by default", {
  # process.autocorr()'s default max_lag is the shortest raw trial length
  # across the list (in contrast to compute.autocorr()'s default, which uses
  # each neuron's own trial length) -- this is why row16$tau above differs
  # very slightly from the standalone compute.autocorr() value pinned above.
  neurons <- dg_demo_neurons()
  trial_lengths <- sapply(neurons, function(nrn) nrow(nrn$fetch_trial_data_R()))
  expect_equal(min(trial_lengths), 1420)
})

test_that("compute.autocorr() default max_lag uses the neuron's own trial length", {
  nrn <- new.neuron(t_per_bin = 10.0)
  raster <- make_synthetic_raster()
  raster <- raster[raster$cell == 1, ]
  neurons <- load.rasters.as.neurons(raster, bin_size = 10.0)
  n1 <- neurons[["neuron_1"]]

  own_length <- nrow(n1$fetch_trial_data_R())
  compute.autocorr(n1, use_raw = TRUE)
  expect_length(n1$fetch_autocorr_R(), own_length %/% 10)
})
