# Regression values captured by running the "Visualizing cross-correlation"
# section of vignettes/tutorial_synchrony_crosscorr.Rmd against
# inst/extdata/cross_corr_demo_data.csv with set.seed(12345).

test_that("compute.crosscorr() matches pinned values for the vignette demo pair", {
  raster <- crosscorr_demo_raster()
  trial_duration <- 100; lag_divider <- 2; bs <- 1.0; pre_ratio <- 0.25
  ref_cell <- 13; comp_cell <- 15

  raster_demo <- make_trials_from_spikes(
    raster = raster,
    cell_num = ref_cell,
    prespike_time = trial_duration * pre_ratio,
    postspike_time = trial_duration * (1 - pre_ratio)
  )
  expect_equal(length(unique(raster_demo$trial)), 156)

  neurons <- load.rasters.as.neurons(
    raster_demo,
    bin_size = bs,
    min_duration = trial_duration,
    max_displacement = -trial_duration * pre_ratio
  )
  ref_neuron <- neurons[[paste0("neuron_", ref_cell)]]
  comp_neuron <- neurons[[paste0("neuron_", comp_cell)]]

  cross_corr <- compute.crosscorr(
    ref_neuron, comp_neuron, "sum", trial_duration / lag_divider, TRUE, FALSE
  )

  expect_length(cross_corr, 50)
  expect_equal(
    head(cross_corr, 5),
    c(0.0001909855, 0.0001273237, 0.0001273237, 0.0004456328, 0.0001273237),
    tolerance = 1e-6
  )
})

test_that("loading lagged spikes produces a cross-correlation peak at the loaded lag", {
  # vignettes/tutorial_synchrony_crosscorr.Rmd loads 10ms-lagged spikes into
  # the comparison neuron and claims "there is a clear peak ... at a lag of
  # 10ms". Index i in the returned vector corresponds to lag (i - 1) ms here
  # (bin size 1ms), so a peak at lag 10ms is index 11.
  withr_seed <- function(seed, code) {
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    old <- if (has_seed) get(".Random.seed", envir = globalenv()) else NULL
    on.exit({
      if (has_seed) {
        assign(".Random.seed", old, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    })
    set.seed(seed)
    force(code)
  }

  raster <- crosscorr_demo_raster()
  trial_duration <- 100; lag_divider <- 2; bs <- 1.0; pre_ratio <- 0.25
  ref_cell <- 13; comp_cell <- 15

  result <- withr_seed(12345, {
    raster_demo <- make_trials_from_spikes(
      raster = raster, cell_num = ref_cell,
      prespike_time = trial_duration * pre_ratio,
      postspike_time = trial_duration * (1 - pre_ratio)
    )

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
    compute.crosscorr(
      ref_neuron_loaded, comp_neuron_loaded, "sum", trial_duration / lag_divider, TRUE, FALSE
    )
  })

  expect_equal(which.max(result), 11)
})

test_that("compute.crosscorr() of a neuron against a copy of itself matches pinned values", {
  # Note: this is close to, but NOT identical to, that neuron's own
  # autocorrelation (compute_autocorrelation uses a different normalization /
  # valid-lag-step convention) -- so this pins compute.crosscorr()'s
  # own output directly rather than asserting equality between the two.
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  nrn_copy <- new.neuron(t_per_bin = fetch.id.data(nrn)$t_per_bin)
  load.trial.data(nrn_copy, fetch.trial.data(nrn))

  cross <- compute.crosscorr(nrn, nrn_copy, "sum", 100, TRUE, FALSE)

  expect_length(cross, 10)
  expect_equal(
    cross,
    c(0.1738209, 0.1258781, 0.1020540, 0.0891096, 0.0808007,
      0.0744066, 0.0692147, 0.0651252, 0.0613976, 0.0578283),
    tolerance = 1e-4
  )
})
