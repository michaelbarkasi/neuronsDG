# S3 methods added to fix CRAN's "S3 generic/method consistency" NOTE (plot.raster/
# plot.autocorrelation looked like unregistered S3 methods for classes "raster"/
# "autocorrelation", but actually took a neuron object). neuron objects are Rcpp Modules
# reference objects -- formally S4 (class "Rcpp_neuron") -- so methods target that class
# directly rather than a wrapper class (prepending a class via class(x) <- c(...) strips
# their S4-ness and breaks the $ accessor Rcpp modules rely on; see git history). The three
# list-returning analysis functions are likewise tagged with their own (plain S3) classes.

test_that("new.neuron() and load.rasters.as.neurons() produce dispatchable neuron objects", {
  nrn <- new.neuron(id_num = 1)
  expect_true(inherits(nrn, "Rcpp_neuron"))
  expect_true(isS4(nrn)) # must stay a formal S4 object for $ method dispatch to keep working

  raster_df <- make_synthetic_raster()
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)
  expect_true(all(vapply(neurons, inherits, logical(1), what = "Rcpp_neuron")))
  expect_true(all(vapply(neurons, isS4, logical(1))))
})

test_that("auto-printing a neuron at the console (show()) delegates to print.Rcpp_neuron()", {
  nrn <- new.neuron(id_num = 9, recording_name = "rec9")
  out <- capture.output(methods::show(nrn))
  expect_match(paste(out, collapse = "\n"), "#9")
  expect_match(paste(out, collapse = "\n"), "rec9")
})

test_that("print.neuron() reports id metadata and reflects what has been computed", {
  nrn <- new.neuron(id_num = 7, recording_name = "rec7", genotype = "KO")
  out_before <- capture.output(print(nrn))
  expect_match(paste(out_before, collapse = "\n"), "#7")
  expect_match(paste(out_before, collapse = "\n"), "rec7")
  expect_match(paste(out_before, collapse = "\n"), "KO")
  expect_match(paste(out_before, collapse = "\n"), "trial data loaded: FALSE")

  neurons <- dg_demo_neurons()
  loaded_nrn <- neurons[[16]]
  out_after <- capture.output(print(loaded_nrn))
  expect_match(paste(out_after, collapse = "\n"), "trial data loaded: TRUE")
  expect_match(paste(out_after, collapse = "\n"), "lambda")

  invisible(capture.output(ret <- print(nrn)))
  expect_identical(ret, nrn)
})

test_that("summary.neuron() reflects whether autocorrelation/fit have been computed", {
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]

  summ_before <- summary(nrn)
  expect_s3_class(summ_before, "summary.neuron")
  expect_false(summ_before$has_autocorr)
  expect_false(summ_before$has_fit)

  compute.autocorr(nrn)
  fit.edf.autocorr(nrn)
  summ_after <- summary(nrn)
  expect_true(summ_after$has_autocorr)
  expect_true(summ_after$has_fit)
  expect_true(all(c("A", "tau", "bias_term") %in% names(summ_after$edf_parameters)))

  out <- capture.output(print(summ_after))
  expect_match(paste(out, collapse = "\n"), "tau")
})

test_that("print.neuronsDG_estimates() summarizes an estimate.autocorr.params() output", {
  neurons <- dg_demo_neurons()[1:2]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE, seed = 1
  )
  expect_s3_class(ests, "neuronsDG_estimates")

  out <- capture.output(print(ests))
  expect_match(paste(out, collapse = "\n"), "neuronsDG_estimates")
  expect_match(paste(out, collapse = "\n"), "mean tau")
})

test_that("print.neuronsDG_analysis() and plot.neuronsDG_analysis() work on an analyze.autocorr() output", {
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE, seed = 1
  )
  result <- analyze.autocorr(ests, covariate = "hemi", n_bs = 20)
  expect_s3_class(result, "neuronsDG_analysis")

  out <- capture.output(print(result))
  expect_match(paste(out, collapse = "\n"), "neuronsDG_analysis")

  expect_identical(plot(result), result)
})

test_that("process.autocorr() output is classed, and plot()/test.sigma.assumption() both work on it", {
  neurons <- dg_demo_neurons()[1:3]
  batch <- process.autocorr(neurons)
  expect_s3_class(batch, "neuronsDG_autocorr_table")
  expect_s3_class(batch, "data.frame")

  expect_no_error(plot(batch))
  expect_no_error(test.sigma.assumption(batch))
  # test.sigma.assumption() also works on a plain, unclassed data frame
  plain_batch <- batch
  class(plain_batch) <- "data.frame"
  expect_no_error(test.sigma.assumption(plain_batch))
})
