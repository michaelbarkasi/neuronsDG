# Replay of vignettes/tutorial_synchrony_crosscorr.Rmd. The single-pair demo
# ("Visualizing cross-correlation") is fast and replayed in full here and in
# test-crosscorrelation.R. The "Analyzing covariates" section's
# cross_corr_by_genotype(), however, loops pairwise cross-correlation over
# every cell in a 520k-row raster for 3 different trial durations and takes
# several minutes PER call (observed: >3 min for just the fastest, 20ms,
# case) -- too slow to replay as a routine test. Its building blocks
# (make_trials_from_spikes() equivalent via load.rasters.as.neurons(), and
# compute_crosscorrelation_R()) are covered directly in test-crosscorrelation.R
# instead.

test_that("the synchrony vignette's single-pair demo matches its own narrative claims", {
  set.seed(12345)
  raster <- crosscorr_demo_raster()

  trial_duration <- 100; lag_divider <- 2; bs <- 1.0; pre_ratio <- 0.25
  ref_cell <- 13; comp_cell <- 15

  raster_demo <- make_trials_from_spikes(
    raster = raster, cell_num = ref_cell,
    prespike_time = trial_duration * pre_ratio,
    postspike_time = trial_duration * (1 - pre_ratio)
  )
  neurons <- load.rasters.as.neurons(
    raster_demo, bin_size = bs,
    min_duration = trial_duration, max_displacement = -trial_duration * pre_ratio
  )
  ref_neuron <- neurons[[paste0("neuron_", ref_cell)]]
  comp_neuron <- neurons[[paste0("neuron_", comp_cell)]]

  cross_corr <- ref_neuron$compute_crosscorrelation_R(
    comp_neuron, "sum", trial_duration / lag_divider, TRUE, FALSE
  )
  # "these two neurons seem to have no obvious peak in their
  #  cross-correlation" -- i.e. a fairly flat profile, not a sharp spike
  expect_lt(max(cross_corr) / mean(cross_corr), 5)

  load_lagged_spikes <- function(raster, cell_num_comp, delay_time) {
    cell_idx_comp <- which(raster$cell == cell_num_comp)
    cell_idx_comp <- sample(cell_idx_comp, length(cell_idx_comp) / 4, replace = FALSE)
    raster$time_in_ms[cell_idx_comp] <- delay_time + rnorm(length(cell_idx_comp), 0, 1)
    raster
  }
  raster_demo_loaded <- load_lagged_spikes(raster_demo, comp_cell, 10)
  neurons_loaded <- load.rasters.as.neurons(
    raster_demo_loaded, bin_size = bs,
    min_duration = trial_duration, max_displacement = -trial_duration * pre_ratio
  )
  ref_neuron_loaded <- neurons_loaded[[paste0("neuron_", ref_cell)]]
  comp_neuron_loaded <- neurons_loaded[[paste0("neuron_", comp_cell)]]
  cross_corr_loaded <- ref_neuron_loaded$compute_crosscorrelation_R(
    comp_neuron_loaded, "sum", trial_duration / lag_divider, TRUE, FALSE
  )

  # "Notice that now there is a clear peak ... at a lag of 10ms"
  expect_equal(which.max(cross_corr_loaded), 11) # index 11 == lag 10ms (1ms bins)
  expect_gt(max(cross_corr_loaded) / mean(cross_corr_loaded), 5)

  # "The above histogram shows a clear peak at 10ms"
  loaded_raster <- comp_neuron_loaded$fetch_spike_raster_R()
  loaded_spikes <- loaded_raster[, 1]
  hist_counts <- hist(loaded_spikes, breaks = 50, plot = FALSE)
  peak_mid <- hist_counts$mids[which.max(hist_counts$counts)]
  expect_true(abs(peak_mid - 10) < diff(hist_counts$breaks)[1])
})

test_that("cross_corr_by_genotype() over the full demo data is not replayed (too slow for routine testing)", {
  skip(paste(
    "cross_corr_by_genotype() takes multiple minutes per call on the full",
    "~23-cell, 520k-row demo raster; its building blocks are covered",
    "directly in test-crosscorrelation.R and the single-pair demo above."
  ))
})
