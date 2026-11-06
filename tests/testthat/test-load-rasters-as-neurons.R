test_that("load.rasters.as.neurons() creates one neuron per unique cell", {
  raster_df <- make_synthetic_raster()
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)

  expect_length(neurons, 3)
  expect_equal(names(neurons), c("neuron_1", "neuron_2", "neuron_3"))
})

test_that("load.rasters.as.neurons() propagates metadata columns", {
  raster_df <- make_synthetic_raster()
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)
  id <- fetch.id.data(neurons[["neuron_1"]])

  expect_equal(id$id_num, 1)
  expect_equal(id$recording_name, "synthetic_rec")
  expect_equal(id$genotype, "WT")
  expect_equal(id$sex, "female")
  expect_equal(id$hemi, "left")
  expect_equal(id$region, "V1")
  expect_equal(id$age, "adult")
  expect_equal(id$unit_time, "ms")
  expect_equal(id$unit_data, "spike")
  expect_equal(id$t_per_bin, 10.0)
})

test_that("load.rasters.as.neurons() falls back to 'not_provided' metadata when absent", {
  raster_df <- make_synthetic_raster()[, c("cell", "time_in_ms", "trial")]
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)
  id <- fetch.id.data(neurons[["neuron_1"]])

  expect_equal(id$recording_name, "not_provided")
  expect_equal(id$genotype, "not_provided")
  expect_equal(id$sex, "not_provided")
  expect_equal(id$hemi, "not_provided")
  expect_equal(id$region, "not_provided")
  expect_equal(id$age, "not_provided")
})

test_that("load.rasters.as.neurons() builds the spike raster each neuron reports", {
  raster_df <- make_synthetic_raster()
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)

  sr1 <- fetch.spike.raster(neurons[["neuron_1"]])
  expect_equal(nrow(sr1), 8) # 4 spikes in trial 1 + 4 in trial 2
  expect_setequal(sr1[, "time"], c(10, 20, 30, 40, 15, 25, 35, 45))
  expect_setequal(sr1[, "trial"], c(1, 2))

  sr2 <- fetch.spike.raster(neurons[["neuron_2"]])
  expect_equal(nrow(sr2), 4)

  sr3 <- fetch.spike.raster(neurons[["neuron_3"]])
  expect_equal(nrow(sr3), 2)
})

test_that("load.rasters.as.neurons() respects time_cutoff", {
  # Restricted to cells 1-2: cell 3 has only one spike surviving a 30ms
  # cutoff (-100 and 100ms), and a single remaining spike cannot define a
  # trial duration (a separate, pre-existing edge case, not under test here).
  raster_df <- make_synthetic_raster()
  raster_df <- raster_df[raster_df$cell %in% c(1, 2), ]
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0, time_cutoff = 30)

  sr1 <- fetch.spike.raster(neurons[["neuron_1"]])
  # Spikes at 10, 20, 30 (trial 1) and 15, 25 (trial 2) survive a 30ms cutoff;
  # 40 and 35, 45 do not.
  expect_true(all(sr1[, "time"] <= 30))
  expect_equal(nrow(sr1), 5)
})

test_that("a single spike left after time_cutoff raises an error", {
  # Documents existing behavior: filtering down to one spike for a cell means
  # no trial duration can be inferred, so load.rasters.as.neurons() errors.
  raster_df <- make_synthetic_raster()
  raster_df <- raster_df[raster_df$cell == 3, ]
  expect_error(
    load.rasters.as.neurons(raster_df, bin_size = 10.0, time_cutoff = 30),
    "zero or negative duration"
  )
})

test_that("load.rasters.as.neurons() accepts a csv file path, not just a data.frame", {
  raster_df <- make_synthetic_raster()
  tmp <- tempfile(fileext = ".csv")
  write.csv(raster_df, tmp, row.names = FALSE)
  on.exit(unlink(tmp))

  neurons_from_path <- load.rasters.as.neurons(tmp, bin_size = 10.0)
  neurons_from_df   <- load.rasters.as.neurons(raster_df, bin_size = 10.0)

  expect_equal(names(neurons_from_path), names(neurons_from_df))
  expect_equal(
    fetch.spike.raster(neurons_from_path[["neuron_1"]]),
    fetch.spike.raster(neurons_from_df[["neuron_1"]])
  )
})

test_that("load.rasters.as.neurons() on real demo data matches pinned structure (regression)", {
  neurons <- dg_demo_neurons()

  expect_length(neurons, 41)
  expect_equal(names(neurons)[c(1, 16)], c("neuron_1", "neuron_16"))
  expect_equal(nrow(fetch.trial.data(neurons[[16]])), 1540)
  expect_equal(ncol(fetch.trial.data(neurons[[16]])), 323)
})
