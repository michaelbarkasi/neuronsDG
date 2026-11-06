# Regression values captured with seed = 42, a 3-neuron subset of the demo spike-raster data,
# n_trials_per_sim = 50, n_sims_per_neurons = 5. Kept small so this file runs in seconds; see
# test-vignette-tau-est-DG.R and test-vignette-tau-est-kilosort.R for the full-fidelity (slow)
# replays using the tutorials' actual parameters.
#
# estimate.autocorr.params() previously crashed on exactly this kind of call
# (realistic multi-neuron data, default bin size) with "max_lag must be less
# than or equal to the trial length" -- a units-conversion bug across
# R/neuronsDG.R's max_lag default and src/neuron.cpp's
# neuron::estimate_autocorr_params(). These tests guard against that
# regression as much as they check the numeric output.
#
# estimate.autocorr.params()'s simulations are drawn by their own C++ random number generator
# (see src/neuron.cpp's dg_mvn_draw_cholesky()), not R's -- set.seed() has no effect on them,
# only this function's own seed argument does. analyze.autocorr()'s bootstrap, in contrast, is
# plain R sample(), so it's still controlled by set.seed() as before.

test_that("estimate.autocorr.params()'s simulated lambda_bin tracks the neuron's true firing rate", {
  # Regression guard for a bug found while investigating how different the native-RNG port's
  # results were from the old implementation's (neither is "noise" -- see git history): inside
  # dg_simulation(), load_trial_data() (which computes lambda_bin = lambda * t_per_bin as it
  # loads) was being called *before* t_per_bin was reset from this neuron's own value (e.g. 10)
  # to the simulated neuron's bin-unit value (1), inflating simulated lambda_bin by that factor.
  # tau/A were unaffected (they come from lambda and the trial-data matrix directly, not
  # lambda_bin), but the lambda_bin/lambda_ms columns estimate.autocorr.params() reports were
  # wrong for any neuron whose own t_per_bin != 1 -- i.e. essentially always, since 1 isn't the
  # default.
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  true_lambda_bin <- unname(fetch.lambda(nrn)[2])
  expect_gt(true_lambda_bin, 0.01) # sanity: this neuron's t_per_bin (10) actually matters here

  ests <- estimate.autocorr.params(
    neuron_list = list(nrn), n_trials_per_sim = 300, n_sims_per_neurons = 30,
    use_raw = TRUE, seed = 7
  )
  expect_equal(mean(ests$estimates$lambda_bin), true_lambda_bin, tolerance = 0.1)
})

test_that("estimate.autocorr.params() runs without error on a realistic multi-neuron list", {
  neurons <- dg_demo_neurons()[1:3]

  expect_no_error(
    ests <- estimate.autocorr.params(
      neuron_list = neurons,
      n_trials_per_sim = 50,
      n_sims_per_neurons = 5,
      use_raw = TRUE,
      seed = 42
    )
  )
  expect_equal(dim(ests$estimates), c(15, 9))
  expect_equal(nrow(ests$neuron_id), 3)
  expect_equal(ests$n_sims_per_neurons, 5)
})

test_that("estimate.autocorr.params() matches pinned regression values", {
  # estimate.autocorr.params() fits an EDF model (via NLopt) to simulated autocorrelation
  # curves derived from each real neuron's own parameters; that fit is not bit-reproducible
  # across compilers/math libraries for borderline/noisy cases (known, accepted behavior --
  # see neuron::fit_autocorrelation() in src/neuron.cpp), so skip on CRAN's diverse build
  # machines rather than pin exact values that can legitimately differ by platform.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons,
    n_trials_per_sim = 50,
    n_sims_per_neurons = 5,
    use_raw = TRUE,
    seed = 42
  )

  expect_equal(
    colnames(ests$estimates),
    c("lambda_ms", "lambda_bin", "A", "tau", "bias_term",
      "autocorr1", "max_autocorr", "mean_autocorr", "min_autocorr")
  )

  # First simulation of each of the 3 neurons (rows 1, 6, 11)
  expect_equal(ests$estimates$A[1],   0.0172710816, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[1], 41.6892651244, tolerance = 1e-2)
  expect_equal(ests$estimates$A[6],   0.0040958268, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[6], 37.0979091081, tolerance = 1e-2)
  expect_equal(ests$estimates$A[11],  0.0191032390, tolerance = 1e-5)
  expect_equal(ests$estimates$tau[11], 58.4879648201, tolerance = 1e-2)
})

test_that("estimate.autocorr.params() is deterministic given the same seed, and differs otherwise", {
  neurons <- dg_demo_neurons()[1:3]
  ests1 <- estimate.autocorr.params(neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, seed = 5)
  ests2 <- estimate.autocorr.params(neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, seed = 5)
  ests3 <- estimate.autocorr.params(neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, seed = 6)
  expect_equal(ests1$estimates, ests2$estimates)
  expect_false(isTRUE(all.equal(ests1$estimates, ests3$estimates)))
})

test_that("estimate.autocorr.params() is unaffected by R's set.seed()", {
  neurons <- dg_demo_neurons()[1:3]
  set.seed(1)
  ests1 <- estimate.autocorr.params(neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, seed = 5)
  set.seed(99999)
  ests2 <- estimate.autocorr.params(neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, seed = 5)
  expect_equal(ests1$estimates, ests2$estimates)
})

test_that("estimate.autocorr.params() works across neurons with heterogeneous trial lengths", {
  # The max_lag bug this guards against was masked for a single neuron and only
  # surfaced with >=2 neurons of differing raw trial length -- which is the
  # realistic case (see test-autocorrelation.R's min(trial_lengths) check).
  neurons <- dg_demo_neurons()
  trial_lengths <- sapply(neurons, function(nrn) nrow(fetch.trial.data(nrn)))
  expect_gt(length(unique(trial_lengths)), 1) # sanity: data really is heterogeneous

  shortest_and_longest <- neurons[c(which.min(trial_lengths), which.max(trial_lengths))]
  expect_no_error(
    estimate.autocorr.params(
      neuron_list = shortest_and_longest,
      n_trials_per_sim = 20,
      n_sims_per_neurons = 2,
      use_raw = TRUE,
      seed = 1
    )
  )
})

test_that("analyze.autocorr() bootstraps tau by covariate and matches pinned values", {
  # See "estimate.autocorr.params() matches pinned regression values" above for why the
  # pinned values below are platform-sensitive and CRAN-skipped.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 50, n_sims_per_neurons = 5, use_raw = TRUE, seed = 42
  )

  set.seed(42) # controls analyze.autocorr()'s own bootstrap sample(); unrelated to the seed above
  result <- analyze.autocorr(ests, covariate = "hemi", n_bs = 100)

  expect_named(result, c("resamples", "distribution_plot"))
  expect_equal(colnames(result$resamples), c("LH", "RH"))
  expect_equal(nrow(result$resamples), 100)
  expect_equal(result$resamples$LH[1], 41.3102831924, tolerance = 1e-2)
  expect_equal(result$resamples$RH[1], 43.7148102970, tolerance = 1e-2)
  expect_s3_class(result$distribution_plot, "ggplot")
})

test_that("analyze.autocorr() rejects input that isn't an estimate.autocorr.params() output", {
  expect_error(analyze.autocorr(list(foo = 1), covariate = "hemi"), "not a valid")
  expect_error(analyze.autocorr(data.frame(x = 1), covariate = "hemi"), "not a valid")
})

test_that("analyze.autocorr() rejects an unknown covariate", {
  # The covariate-name check itself doesn't depend on real data, but building
  # `ests` does: estimate.autocorr.params() fits each neuron's own EDF model
  # then simulates from it via a Cholesky-based multivariate draw (see
  # dg_mvn_draw_cholesky() in src/neuron.cpp), and a poorly-converged fit on
  # real, noisy data (the same platform-sensitive NLopt fragility noted above)
  # can feed that draw a near-singular matrix and throw "matrix seems negative
  # semi-definite" before this test ever reaches the check it's meant to
  # exercise. Skip on CRAN for the same reason as the pinned-value tests.
  skip_on_cran()
  skip_on_ci()
  neurons <- dg_demo_neurons()[1:3]
  ests <- estimate.autocorr.params(
    neuron_list = neurons, n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE, seed = 42
  )
  expect_error(analyze.autocorr(ests, covariate = "not_a_real_covariate"), "not found")
})

test_that("analyze.autocorr() accepts a list of estimate.autocorr.params() outputs (batched input)", {
  neurons <- dg_demo_neurons()[1:3]
  ests1 <- estimate.autocorr.params(
    neuron_list = neurons[1:2], n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE, seed = 42
  )
  ests2 <- estimate.autocorr.params(
    neuron_list = neurons[3], n_trials_per_sim = 20, n_sims_per_neurons = 2, use_raw = TRUE, seed = 43
  )

  expect_no_error(
    result <- analyze.autocorr(list(ests1, ests2), covariate = "hemi", n_bs = 20)
  )
  expect_true(is.data.frame(result$resamples))
})
