# Regression values captured with set.seed(42), a 3-neuron subset of
# inst/extdata/spike_rasters_2022data.csv, n_trials_per_sim = 50,
# n_sims_per_neurons = 5. Kept small so this file runs in seconds; see
# test-vignette-tau-est-DG.R and test-vignette-tau-est-kilosort.R for the
# full-fidelity (slow) replays using the tutorials' actual parameters.
#
# estimate.autocorr.params() previously crashed on exactly this kind of call
# (realistic multi-neuron data, default bin size) with "max_lag must be less
# than or equal to the trial length" -- a units-conversion bug across
# R/neuronsDG.R's max_lag default and src/neuron.cpp's
# neuron::estimate_autocorr_params(). These tests guard against that
# regression as much as they check the numeric output.

test_that("estimate.autocorr.params() runs without error on a realistic multi-neuron list", {
  set.seed(42)
  neurons <- dg_demo_neurons()[1:3]

  expect_no_error(
    ests <- estimate.autocorr.params(
      neuron_list = neurons,
      n_trials_per_sim = 50,
      n_sims_per_neurons = 5,
      use_raw = TRUE
    )
  )
  expect_equal(dim(ests$estimates), c(15, 9))
  expect_equal(nrow(ests$neuron_id), 3)
  expect_equal(ests$n_sims_per_neurons, 5)
})

test_that("estimate.autocorr.params() matches pinned regression values", {
  set.seed(42)
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons,
    n_trials_per_sim = 50,
    n_sims_per_neurons = 5,
    use_raw = TRUE
  )

  expect_equal(
    colnames(ests$estimates),
    c("lambda_ms", "lambda_bin", "A", "tau", "bias_term",
      "autocorr1", "max_autocorr", "mean_autocorr", "min_autocorr")
  )

  # First simulation of each of the 3 neurons (rows 1, 6, 11)
  expect_equal(ests$estimates$A[1],   0.0084601196, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[1], 24.2421343255, tolerance = 1e-2)
  expect_equal(ests$estimates$A[6],   0.0314380844, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[6], 50.2193586938, tolerance = 1e-2)
  expect_equal(ests$estimates$A[11],  0.0310142998, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[11], 47.1933049604, tolerance = 1e-2)
})

test_that("estimate.autocorr.params() works across neurons with heterogeneous trial lengths", {
  # The bug this guards against was masked for a single neuron and only
  # surfaced with >=2 neurons of differing raw trial length -- which is the
  # realistic case (see test-autocorrelation.R's min(trial_lengths) check).
  set.seed(1)
  neurons <- dg_demo_neurons()
  trial_lengths <- sapply(neurons, function(nrn) nrow(nrn$fetch_trial_data_R()))
  expect_gt(length(unique(trial_lengths)), 1) # sanity: data really is heterogeneous

  shortest_and_longest <- neurons[c(which.min(trial_lengths), which.max(trial_lengths))]
  expect_no_error(
    estimate.autocorr.params(
      neuron_list = shortest_and_longest,
      n_trials_per_sim = 20,
      n_sims_per_neurons = 2,
      use_raw = TRUE
    )
  )
})

test_that("analyze.autocorr() bootstraps tau by covariate and matches pinned values", {
  set.seed(42)
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 50, n_sims_per_neurons = 5, use_raw = TRUE
  )

  set.seed(42)
  result <- analyze.autocorr(ests, covariate = "hemi", n_bs = 100)

  expect_named(result, c("resamples", "distribution_plot"))
  expect_equal(colnames(result$resamples), c("LH", "RH"))
  expect_equal(nrow(result$resamples), 100)
  expect_equal(result$resamples$LH[1], 27.94411, tolerance = 1e-2)
  expect_equal(result$resamples$RH[1], 49.22139, tolerance = 1e-2)
  expect_s3_class(result$distribution_plot, "ggplot")
})

test_that("analyze.autocorr() rejects input that isn't an estimate.autocorr.params() output", {
  expect_error(analyze.autocorr(list(foo = 1), covariate = "hemi"), "not a valid")
  expect_error(analyze.autocorr(data.frame(x = 1), covariate = "hemi"), "not a valid")
})

test_that("analyze.autocorr() rejects an unknown covariate", {
  set.seed(42)
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE
  )
  expect_error(analyze.autocorr(ests, covariate = "not_a_real_covariate"), "not found")
})

test_that("analyze.autocorr() accepts a list of estimate.autocorr.params() outputs (batched input)", {
  set.seed(42)
  neurons <- dg_demo_neurons()[1:3]
  ests1 <- estimate.autocorr.params(
    neuron_list = neurons[1:2], n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE
  )
  ests2 <- estimate.autocorr.params(
    neuron_list = neurons[3], n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE
  )

  expect_no_error(
    result <- analyze.autocorr(list(ests1, ests2), covariate = "hemi", n_bs = 20)
  )
  expect_true(is.data.frame(result$resamples))
})
