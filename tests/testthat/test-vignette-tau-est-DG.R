# Full-fidelity replay of vignettes/articles/tutorial_tau_est_DG.Rmd: same data, same seed,
# same parameters as the tutorial itself (not a scaled-down version). Slow (~55s dominated by
# estimate.autocorr.params()'s 4,100 DG simulations), so skipped on CRAN; run locally / in CI
# to catch regressions in either the numbers or the narrative claims the tutorial makes about
# them.
#
# estimate.autocorr.params()'s simulations are drawn by their own C++ random number generator,
# not R's -- set.seed() below has no effect on them, only the seed argument passed to that call
# does. analyze.autocorr()'s bootstrap, in contrast, is plain R sample(), so it's still
# controlled by set.seed() as before.

test_that("the DG tutorial's full pipeline reproduces its own narrative claims and pinned values", {
  skip_on_cran()
  skip_on_ci() # ~80s; keep CI fast, run this one locally / explicitly instead

  set.seed(12345)
  spike.rasters <- get_spike_rasters_2022_data()
  if (is.null(spike.rasters)) skip("spike_rasters_2022data.csv demo data unavailable (no internet?)")
  neurons <- load.rasters.as.neurons(spike.rasters, bin_size = 10.0)
  expect_length(neurons, 41)

  cell_high <- 16; cell_low <- 1

  # "The first has high autocorrelation ... The second has low autocorrelation"
  raster_high <- as.data.frame(fetch.spike.raster(neurons[[cell_high]]))
  raster_low  <- as.data.frame(fetch.spike.raster(neurons[[cell_low]]))

  compute.autocorr(neurons[[cell_high]], use_raw = TRUE)
  fit.edf.autocorr(neurons[[cell_high]])
  compute.autocorr(neurons[[cell_low]], use_raw = TRUE)
  fit.edf.autocorr(neurons[[cell_low]])

  edf_high <- fetch.edf.parameters(neurons[[cell_high]])
  edf_low  <- fetch.edf.parameters(neurons[[cell_low]])

  expect_equal(unname(edf_high["tau"]), 69.3777201221, tolerance = 1e-3)
  expect_equal(unname(edf_low["tau"]),  37.4334419605, tolerance = 1e-3)
  # "the time constant tau is estimated to be Xms and the initial
  #  autocorrelation A is estimated to be Y" (narrative, high-autocorr cell)
  expect_equal(round(unname(edf_high["tau"]), 1), 69.4)
  expect_equal(round(unname(edf_high["A"]), 3), 0.127)

  batch <- process.autocorr(neurons)
  expect_equal(nrow(batch), 41)

  # "For the purpose of speed, this tutorial runs only 100 simulations per
  # neuron" / n_trials_per_sim = 500
  autocor.ests <- estimate.autocorr.params(
    neuron_list = neurons,
    n_trials_per_sim = 500,
    n_sims_per_neurons = 100,
    use_raw = TRUE,
    seed = 12345
  )
  expect_equal(dim(autocor.ests$estimates), c(4100, 9))
  expect_equal(autocor.ests$estimates$A[1],   0.0116827377, tolerance = 1e-4)
  expect_equal(autocor.ests$estimates$tau[1], 37.6241013621, tolerance = 1e-1)

  autocor.results.bootstraps <- analyze.autocorr(
    autocor.ests,
    covariate = c("hemi", "genotype"),
    n_bs = 1e4
  )
  resamples <- autocor.results.bootstraps$resamples
  expect_setequal(colnames(resamples), c("LH_CBA/J", "RH_CBA/J"))

  # "reproduces [Neophytou et al. 2022's] analysis" that the right auditory
  # cortex has MORE recurrence (higher tau) than the left
  means <- colMeans(resamples)
  expect_gt(means[["RH_CBA/J"]], means[["LH_CBA/J"]])
  expect_equal(means[["LH_CBA/J"]], 54.196845, tolerance = 1)
  expect_equal(means[["RH_CBA/J"]], 57.963491, tolerance = 1)

  expect_s3_class(autocor.results.bootstraps$distribution_plot, "ggplot")
})
