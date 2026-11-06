# Full-fidelity replay of vignettes/articles/tutorial_tau_est_kilosort.Rmd: same real kilosort4
# demo data, same seed, same parameters as the tutorial (not a scaled-down version). Slow
# (dominated by estimate.autocorr.params()'s 3,600 DG simulations, on top of the ~40s
# preprocess.kilo4() shared with test-preprocessing.R via the memoized fixture), so skipped on
# CRAN/CI; run locally / explicitly to catch regressions in either the numbers or the narrative
# claims the tutorial makes about them.
#
# estimate.autocorr.params()'s simulations are drawn by their own C++ random number generator,
# not R's -- set.seed() below has no effect on them, only the seed argument passed to that call
# does. analyze.autocorr()'s bootstrap, in contrast, is plain R sample(), so it's still
# controlled by set.seed() as before.

test_that("the kilosort tutorial's full pipeline reproduces its own narrative claims and pinned values", {
  skip_on_cran()
  skip_on_ci() # slow; run this one locally / explicitly instead

  set.seed(12345)
  spike.rasters <- kilo4_demo_preprocessed()

  neurons <- load.rasters.as.neurons(
    spike.rasters$spikes,
    bin_size = 10.0,
    sample_rt = 1e3
  )
  expect_length(neurons, 37)

  cell_high <- 25; cell_low <- 13

  compute.autocorr(neurons[[cell_high]]); fit.edf.autocorr(neurons[[cell_high]])
  compute.autocorr(neurons[[cell_low]]);  fit.edf.autocorr(neurons[[cell_low]])

  edf_high <- fetch.edf.parameters(neurons[[cell_high]])
  edf_low  <- fetch.edf.parameters(neurons[[cell_low]])

  # "the high-autocorrelation cell has a time constant of Xms and an
  #  amplitude of Y, while the low-autocorrelation cell has a time constant
  #  of only Zms and an amplitude of W"
  expect_gt(unname(edf_high["tau"]), unname(edf_low["tau"]))
  expect_equal(round(unname(edf_high["tau"]), 1), 104.4)
  expect_equal(round(unname(edf_low["tau"]), 1), 44.2)

  batch <- process.autocorr(neurons)
  expect_equal(nrow(batch), 37)

  # "the last neuron, number 37, is the only C57 sample"
  expect_equal(spike.rasters$cluster.key$cell[37], 37)
  last_cov <- spike.rasters$cluster.key[spike.rasters$cluster.key$cell == 37, ]
  expect_equal(last_cov$genotype, "C57")

  neurons <- neurons[-37]
  expect_length(neurons, 36)

  autocor.ests <- estimate.autocorr.params(
    neuron_list = neurons,
    n_trials_per_sim = 500,
    n_sims_per_neurons = 100,
    seed = 12345
  )
  expect_equal(dim(autocor.ests$estimates), c(3600, 9))
  expect_equal(autocor.ests$estimates$A[1],   0.0123871968, tolerance = 1e-4)
  expect_equal(autocor.ests$estimates$tau[1], 59.4834104130, tolerance = 1e-1)

  autocor.results.bootstraps <- analyze.autocorr(
    autocor.ests,
    covariate = c("hemi", "genotype"),
    n_bs = 1e4
  )
  expect_true(is.data.frame(autocor.results.bootstraps$resamples))
  expect_false(anyNA(autocor.results.bootstraps$resamples))
  expect_s3_class(autocor.results.bootstraps$distribution_plot, "ggplot")
})
