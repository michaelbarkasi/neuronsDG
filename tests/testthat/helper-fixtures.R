# Shared fixtures for the neuronsDG test suite.
#
# Heavy fixtures (anything built from the real demo data in inst/extdata) are
# memoized in .fixture_cache so that, however many test files use them, the
# underlying computation (reading a 20MB csv, preprocessing kilosort4 data,
# etc.) only runs once per test session.

.fixture_cache <- new.env(parent = emptyenv())

get_fixture <- function(name, make) {
  if (!exists(name, envir = .fixture_cache, inherits = FALSE)) {
    assign(name, make(), envir = .fixture_cache)
  }
  get(name, envir = .fixture_cache, inherits = FALSE)
}

# A small, fully synthetic compact raster: 3 cells, 4 trials, hand-picked spike
# times so expected results can be worked out by hand rather than just pinned.
make_synthetic_raster <- function() {
  data.frame(
    cell        = c(1,1,1,1, 1,1,1,1, 2,2,2,2, 3,3),
    time_in_ms  = c(10,20,30,40, 15,25,35,45, 12,22,32,42, 100,-100),
    trial       = c(1,1,1,1, 2,2,2,2, 1,1,1,1, 1,1),
    recording_name = "synthetic_rec",
    genotype    = "WT",
    # Deliberately not "F": a single-letter column of all "F" (or all "T")
    # values is silently read back as logical FALSE/TRUE by read.csv()'s type
    # inference, not character -- a base-R csv gotcha, not a package bug, but
    # one the csv-roundtrip test below would otherwise trip over.
    sex         = "female",
    hemi        = "left",
    region      = "V1",
    age         = "adult"
  )
}

# Real demo data used directly by the vignettes; small enough (3.2MB) that
# loading it is fast.
dg_demo_spike_rasters <- function() {
  get_fixture("dg_demo_spike_rasters", function() {
    read.csv(system.file("extdata", "spike_rasters_2022data.csv", package = "neuronsDG"))
  })
}

# neuron list built exactly as in vignettes/tutorial_tau_est_DG.Rmd.
#
# Deliberately NOT memoized: the underlying Rcpp 'neuron' objects are mutable
# reference objects, and functions like compute.autocorr()/fit.edf.autocorr()/
# estimate.autocorr.params() mutate them in place (storing autocorr/autocorr_edf/
# sigma_gauss). Sharing one cached list across tests would let one test's
# mutations leak into another's. The underlying csv read IS memoized above, so
# rebuilding the neuron list here is cheap (<1s).
dg_demo_neurons <- function() {
  load.rasters.as.neurons(dg_demo_spike_rasters(), bin_size = 10.0)
}

# Kilosort4 demo data, preprocessed exactly as in
# vignettes/tutorial_tau_est_kilosort.Rmd. This is the expensive one (~40s):
# reads real kilosort4 output (npy/tsv/mat files) across 7 recording sessions.
kilo4_metadata_formatted <- function() {
  get_fixture("kilo4_metadata_formatted", function() {
    kilo4_metadata <- read.csv(
      system.file("extdata", "meta_data_kilo4demo.csv", package = "neuronsDG")
    )
    rownames(kilo4_metadata) <- paste0(kilo4_metadata$DAY, "_", kilo4_metadata$Neuralynx_ID)
    kilo4_metadata <- kilo4_metadata[, c("HEMISPHERE", "STRAIN", "AGE", "SEX")]
    colnames(kilo4_metadata) <- c("hemi", "genotype", "age", "sex")
    kilo4_metadata
  })
}

kilo4_demo_preprocessed <- function() {
  get_fixture("kilo4_demo_preprocessed", function() {
    preprocess.kilo4(
      trial_time_start = 500,
      trial_time_end = 500 + 1520,
      recording.folder = system.file("extdata", "kilo4demo", package = "neuronsDG"),
      meta_data = kilo4_metadata_formatted(),
      max_spikes = 1e4,
      min_spikes = 1e2,
      min_trials = 1e2,
      pure_trials_only = TRUE,
      good_cells_only = TRUE,
      stim_responsive_only = TRUE,
      verbose = FALSE
    )
  })
}

# neuron list built exactly as in vignettes/tutorial_tau_est_kilosort.Rmd.
# Not memoized, for the same reason as dg_demo_neurons() above; the expensive
# step (preprocess.kilo4(), ~40s) is in kilo4_demo_preprocessed() and stays
# memoized since it returns only plain data frames/matrices.
kilo4_demo_neurons <- function() {
  load.rasters.as.neurons(
    kilo4_demo_preprocessed()$spikes,
    bin_size = 10.0,
    sample_rt = 1e3
  )
}

# Compact raster used by the cross-correlation/synchrony vignette (20MB; only
# read once and shared across tests that need it).
crosscorr_demo_raster <- function() {
  get_fixture("crosscorr_demo_raster", function() {
    read.csv(system.file("extdata", "cross_corr_demo_data.csv", package = "neuronsDG"))
  })
}

# Trial-parsing helper, copied verbatim from
# vignettes/tutorial_synchrony_crosscorr.Rmd, so vignette-replay tests exercise
# exactly the logic shown to users.
make_trials_from_spikes <- function(
    raster, cell_num, prespike_time, postspike_time, clean_only = TRUE
  ) {
  cell_pen <- unique(raster$recording_name[raster$cell == cell_num])
  raster_p <- raster[raster$recording_name == cell_pen, ]
  cell_idx <- which(raster_p$cell == cell_num)
  ref_cell_spikes <- raster_p$time_in_ms[cell_idx]
  ref_cell_spikes_prespike_quiet <- diff(c(0, ref_cell_spikes))
  spikes_preclear <- ref_cell_spikes_prespike_quiet >= prespike_time
  ref_cell_spikes_postspike_quiet <- diff(c(ref_cell_spikes, Inf))
  spikes_postclear <- ref_cell_spikes_postspike_quiet >= postspike_time
  clear_spikes_idx <- cell_idx[spikes_preclear & spikes_postclear]
  raster_p$trial <- 0
  for (i in seq_along(clear_spikes_idx)) {
    spike_time <- raster_p$time_in_ms[clear_spikes_idx[i]]
    trial_start <- spike_time - prespike_time
    trial_end <- spike_time + postspike_time
    trial_mask <- raster_p$time_in_ms >= trial_start & raster_p$time_in_ms < trial_end
    raster_p$trial[trial_mask] <- i
    raster_p$time_in_ms[trial_mask] <- raster_p$time_in_ms[trial_mask] - spike_time
  }
  raster_p <- raster_p[raster_p$trial > 0, ]
  raster_p
}
