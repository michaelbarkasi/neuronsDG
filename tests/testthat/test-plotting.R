test_that("make.plot.title() assembles id info into the title string", {
  nrn <- new.neuron(id_num = 5, recording_name = "rec1", genotype = "KO", sex = "M")
  title <- neuronsDG:::make.plot.title(nrn, "Test title")
  expect_equal(title, "Test title, neuron rec1, #5, KO, M")
})

test_that("make.plot.title() omits covariates left at their 'not_provided' default", {
  nrn <- new.neuron(id_num = 2, recording_name = "rec2")
  title <- neuronsDG:::make.plot.title(nrn, "Raster")
  expect_equal(title, "Raster, neuron rec2, #2")
})

test_that("make.plot.title() flags simulated neurons", {
  nrn <- new.neuron(id_num = 3, recording_name = "rec3", sim = TRUE)
  title <- neuronsDG:::make.plot.title(nrn, "Raster")
  expect_match(title, "#3 \\(sim\\)")
})

test_that("plot.autocorrelation() returns a ggplot built from the neuron's own autocorrelation data", {
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  compute.autocorr(nrn); fit.edf.autocorr(nrn)

  plt <- plot.autocorrelation(nrn, return_plot = TRUE)
  expect_s3_class(plt, "ggplot")
  # Two layers by default: estimated autocorrelation + EDF fit (no bias line
  # unless bias_term is explicitly supplied -- see next test).
  expect_length(plt$layers, 2)

  built <- ggplot2::ggplot_build(plt)
  autocorr <- nrn$fetch_autocorr_R()
  expect_equal(built$data[[1]]$y, autocorr[2:length(autocorr)], tolerance = 1e-8)
})

test_that("plot.autocorrelation()'s documented bias_term default (0) does not match its actual default (NULL)", {
  # man/plot-autocorrelation.Rd documents bias_term's default as 0, but
  # R/neuronsDG.R defines it as NULL. This pins the ACTUAL current default:
  # NULL means no bias-term reference line is drawn.
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  compute.autocorr(nrn); fit.edf.autocorr(nrn)

  plt_default <- plot.autocorrelation(nrn, return_plot = TRUE)
  expect_length(plt_default$layers, 2)

  plt_with_bias <- plot.autocorrelation(nrn, bias_term = 0.05, return_plot = TRUE)
  expect_length(plt_with_bias$layers, 3)
  expect_s3_class(plt_with_bias$layers[[2]]$geom, "GeomHline")
})

test_that("plot.autocorrelation() respects plot_time_cutoff", {
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  compute.autocorr(nrn); fit.edf.autocorr(nrn)

  plt <- plot.autocorrelation(nrn, plot_time_cutoff = 20, return_plot = TRUE)
  built <- ggplot2::ggplot_build(plt)
  expect_true(all(built$data[[1]]$x <= 20))
})

test_that("plot.raster() returns one line segment per spike", {
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]

  plt <- plot.raster(nrn, return_plot = TRUE)
  expect_s3_class(plt, "ggplot")
  expect_length(plt$layers, 1)

  built <- ggplot2::ggplot_build(plt)
  expect_equal(nrow(built$data[[1]]), nrow(nrn$fetch_spike_raster_R()))
})

test_that("plot.raster() draws a stimulus-onset line only when zero_as_onset finds a non-positive time", {
  raster_df <- make_synthetic_raster()
  raster_df <- raster_df[raster_df$cell == 3, ] # includes a spike at -100ms
  neurons <- load.rasters.as.neurons(raster_df, bin_size = 10.0)
  nrn <- neurons[["neuron_3"]]

  plt_with_onset <- plot.raster(nrn, zero_as_onset = TRUE, return_plot = TRUE)
  expect_length(plt_with_onset$layers, 2)

  plt_without_onset <- plot.raster(nrn, zero_as_onset = FALSE, return_plot = TRUE)
  expect_length(plt_without_onset$layers, 1)
})

test_that("plot functions print instead of returning when return_plot = FALSE", {
  neurons <- dg_demo_neurons()
  nrn <- neurons[[16]]
  compute.autocorr(nrn); fit.edf.autocorr(nrn)

  expect_null(suppressWarnings(plot.autocorrelation(nrn, return_plot = FALSE)))
  expect_null(suppressWarnings(plot.raster(nrn, return_plot = FALSE)))
})
