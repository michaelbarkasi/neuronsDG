
# Intro ################################################################################################################

# NeuronsDG: A framework for neuron modelling, simulation, and analysis

# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2026

#' @useDynLib neuronsDG, .registration = TRUE
#' @import Rcpp
NULL

.onLoad <- function(libname, pkgname) {
    Rcpp::loadModule("neuron", TRUE)
    # Our compiled code calls NLopt only through nloptr's registered C routines
    # (nloptrAPI.h, resolved via R_GetCCallable()), not a system libnlopt. Those routines are
    # only registered once nloptr's own namespace has actually been loaded, which merely
    # listing nloptr in Imports does not by itself trigger -- so force it here with a real
    # (harmless) call, which R CMD check also recognizes as actually using the Import.
    invisible(nloptr::is.nloptr(NULL))
    # show()'s S4 method must be registered only after loadModule() (above) has created the
    # "Rcpp_neuron" class -- registering it at the top level of this file instead errors with
    # "no definition for class 'Rcpp_neuron'", since top-level code runs before .onLoad() does.
    # (Dispatches to print.Rcpp_neuron() by name, not generic print(), so R CMD check's static
    # scan of .onLoad() for console-writing calls doesn't flag this.)
    methods::setMethod("show", "Rcpp_neuron", function(object) print.Rcpp_neuron(object))
  }

# Initialization for C++ object classes ################################################################################

#' Initialize neuron
#' 
#' This function initializes a new neuron object with specified parameters. The neuron objects are the main tools for running correlation analyses of spike data.
#' 
#' @param id_num Numeric identifier for the neuron (default: 0).
#' @param recording_name Recording (if any) on which this neuron is based (default: "not_provided").
#' @param type Modeled electrophysiology of neuron, e.g. "generic", "blackbox" "LIF", "McCullochPitts", "excitatory", "inhibitory", etc. (default: "generic").
#' @param genotype Genotype of animal, e.g. "WT", "KO", "MECP2", "transgenic", etc. (default: "not_provided").
#' @param sex Sex of animal (default: "not_provided").
#' @param hemi Hemisphere of neuron in which neuron is located, e.g. "left", "right" (default: "not_provided").
#' @param region Brain region in which neuron is located, e.g. "V1", "M1", "CA1", "PFC", etc. (default: "not_provided").
#' @param age Age of animal, e.g. "P0", "P7", "P14", "adult", etc. (default: "not_provided").
#' @param sim Logical indicating if neuron is simulated (TRUE) or recorded (FALSE) (default: FALSE).
#' @param unit_time Unit of time for spike raster or other recording data, e.g. "ms", "s", etc. (default: "ms").
#' @param unit_sample_rate Unit of sample rate for spike raster or other recording data, e.g. "Hz", "kHz", etc. (default: "Hz").
#' @param unit_data Unit of data for spike raster or other recording data, e.g. "mV", "spike", etc. (default: "mV").
#' @param t_per_bin Time (in above units) per bin, e.g., 1 ms per bin (default: 10.0).
#' @param sample_rate Sample rate (in above units), e.g., 1e4 Hz (default: 1e4).
#' @return A new neuron object.
#' @export
new.neuron <- function(
    id_num = 0, 
    recording_name = "not_provided", 
    type = "generic", 
    genotype = "not_provided",
    sex = "not_provided",
    hemi = "not_provided",
    region = "not_provided",
    age = "not_provided", 
    sim = FALSE, 
    unit_time = "ms", 
    unit_sample_rate = "Hz", 
    unit_data = "mV", 
    t_per_bin = 10.0, 
    sample_rate = 1e4
  ) {
    nrn <- new(
      neuron,
      id_num, recording_name, type, genotype, sex, hemi, region, age, sim, unit_time, unit_sample_rate, unit_data, t_per_bin, sample_rate
    )
    return(nrn)
  }

# Functions for neuron analysis ########################################################################################

#' Load neuron data from raster file
#' 
#' This function loads spike raster data from a data frame or CSV file and converts each unique cell into a neuron object. The raster data frame must contain columns for cell identifier, spike time in milliseconds, and trial number. Optional metadata columns can also be included.
#' 
#' @param raster_df Data frame (or file name to csv importable as such), each row a spike; must have columns: cell, time_in_ms, trial; optional columns: recording_name, hemisphere, genotype, sex, region, age.
#' @param bin_size Size of time bins in milliseconds (default: 10).
#' @param min_duration Minimum trial duration in milliseconds (default: 0).
#' @param max_displacement Max displacement of trial center in milliseconds (default: 1e9, i.e., "infinity").
#' @param sample_rt Sample rate in the default unit for neuron objects, Hz (default: 1e4).
#' @param time_cutoff Maximum time (in ms) to include spikes; spikes occurring after this time will be excluded (default: Inf).
#' @return A list of neuron objects, one per unique cell in the raster data.
#' @export
load.rasters.as.neurons <- function(
    raster_df,
    bin_size = 10.0, 
    min_duration = 0,
    max_displacement = 1e9,
    sample_rt = 1e4,
    time_cutoff = Inf
  ) {
    
    # Load data
    if (is.character(raster_df)) {
      raster_df <- as.data.frame(read.csv(raster_df))
    }
    
    # Initiate list to hold neurons
    neuron_list <- list()
    
    cell_nums <- unique(raster_df$cell)
    for (c in cell_nums) {
      
      # Make name for new neuron
      new_name <- paste0("neuron_", c)
      # Grab raster
      c_mask <- raster_df$cell == c
      raster <- raster_df[c_mask, c("time_in_ms", "trial")]
      colnames(raster) <- c("time", "trial")
      raster <- raster[raster$time <= time_cutoff,]
      # Extract meta data
      recording_name <- "not_provided"
      if ("recording_name" %in% colnames(raster_df)) {
        recording_name <- unique(raster_df$recording_name[c_mask])
        if (length(recording_name) != 1) recording_name <- "not_provided"
      } 
      hemi <- "not_provided"
      if ("hemi" %in% colnames(raster_df)) {
        hemi <- unique(raster_df$hemi[c_mask])
        if (length(hemi) != 1) hemi <- "not_provided"
      }
      genotype <- "not_provided"
      if ("genotype" %in% colnames(raster_df)) {
        genotype <- unique(raster_df$genotype[c_mask])
        if (length(genotype) != 1) genotype <- "not_provided"
      }
      sex <- "not_provided"
      if ("sex" %in% colnames(raster_df)) {
        sex <- unique(raster_df$sex[c_mask])
        if (length(sex) != 1) sex <- "not_provided"
      }
      region <- "not_provided"
      if ("region" %in% colnames(raster_df)) {
        region <- unique(raster_df$region[c_mask])
        if (length(region) != 1) region <- "not_provided"
      }
      age <- "not_provided"
      if ("age" %in% colnames(raster_df)) {
        age <- unique(raster_df$age[c_mask])
        if (length(age) != 1) age <- "not_provided"
      }
      
      # Make new neuron object
      new_neuron <- new.neuron(
        id_num = c,
        recording_name = recording_name,
        genotype = genotype,
        sex = sex,
        hemi = hemi,
        region = region,
        age = age,
        unit_time = "ms",
        unit_data = "spike",
        t_per_bin = bin_size,
        sample_rate = sample_rt
      )
      raster <- as.matrix(raster)
      new_neuron$load_spike_raster_R(raster, min_duration, max_displacement)
      
      # Save in list
      neuron_list[[new_name]] <- new_neuron
      
    }
    
    return(neuron_list)
    
  }

# Helper function for plotting
make.plot.title <- function(
    neuron,                   # neuron object
    plot_title
  ) {
    
    # Grab neuron ID info
    id_data <- neuron$fetch_id_data() 
    recording_name <- id_data$recording_name
    id_num <- id_data$id_num
    genotype <- id_data$genotype
    sex <- id_data$sex
    hemi <- id_data$hemi
    region <- id_data$region
    age <- id_data$age
    sim <- id_data$sim
    if (sim) id_num <- paste(id_num, "(sim)")
    
    # Make plot title 
    plot_title <- paste0(
      plot_title,
      ", neuron ", recording_name, ", #", id_num
    )
    # Add covariates
    for (var in c(genotype, sex, hemi, region, age)) {
      if (var != "not provided" && var != "not_provided") {
        plot_title <- paste0(plot_title, ", ", var)
      }
    }
    
    return(plot_title)
    
  }

#' Plot autocorrelation of a neuron
#'
#' Fetches and plots autocorrelation information from a neuron object. Assumes the autocorrelation has already been computed. Called by \code{plot.neuron(nrn, type = "autocorrelation")}.
#'
#' @param nrn Neuron object for which to plot autocorrelation.
#' @param plot_title Title for the plot (default: "Est. autocorr").
#' @param bias_term Bias term to plot as a horizontal line (default: NULL).
#' @param plot_time_cutoff Maximum lag (in bins) to display on the x-axis (default: Inf).
#' @param return_plot Logical indicating whether to return the ggplot object (TRUE) or print it (FALSE) (default: FALSE).
#' @return A ggplot object showing the estimated and fitted autocorrelation.
#' @export
plot_autocorrelation <- function(
    nrn,
    plot_title = "Est. autocorr",
    bias_term = NULL,
    plot_time_cutoff = Inf,
    return_plot = FALSE
  ) {
    
    # Make title
    plot_title <- make.plot.title(nrn, plot_title)
    
    # Fetch estimated and fitted autocorrelation 
    autocorr <- nrn$fetch_autocorr_R()
    autocorr_edf <- nrn$fetch_autocorr_edf_R()
    
    # Make data frame for plotting 
    df_temp <- data.frame(
      autocorrelation = autocorr[2:length(autocorr)],
      bin = seq(from = 1, to = length(autocorr) - 1, by = 1)
    )
    if (length(autocorr_edf) > 1) {
      df_temp$autocorrelation_fitted <- autocorr_edf[2:length(autocorr_edf)]
    }
    df_temp <- df_temp[df_temp$bin <= plot_time_cutoff,]
    
    # Make plot
    title_size <- 14 
    axis_size <- 12 
    legend_size <- 10
    plt <- ggplot2::ggplot(df_temp) +
      ggplot2::geom_line(ggplot2::aes(x = bin, y = autocorrelation), color = "blue") + 
      ggplot2::labs(
        title = plot_title,
        x = "Lag (bins)",
        y = "Autocorrelation"
      ) + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
        axis.title = ggplot2::element_text(size = axis_size),
        axis.text = ggplot2::element_text(size = axis_size),
        legend.title = ggplot2::element_text(size = legend_size),
        legend.text = ggplot2::element_text(size = legend_size),
        legend.position = "bottom"
      )
    if (!is.null(bias_term)) {
      plt <- plt + ggplot2::geom_hline(yintercept = bias_term, linewidth = 2, linetype = "dotted", color = "darkgray")
    }
    if (length(autocorr_edf) > 1) {
      plt <- plt + ggplot2::geom_line(ggplot2::aes(x = bin, y = autocorrelation_fitted), color = "red")
    }
    if (return_plot) {
      return(plt)
    } else {
      print(plt)
      return(invisible(NULL))
    }
    
  }

#' Plot spike raster of a neuron
#'
#' Fetches and plots the spike raster from a neuron object. Assumes the spike raster has already been loaded. Called by \code{plot.neuron(nrn, type = "raster")} (the default plot for a neuron).
#'
#' @param nrn Neuron object for which to plot spike raster.
#' @param plot_title Title for the plot (default: "Spike raster").
#' @param zero_as_onset Logical indicating whether to plot a vertical line at time zero to indicate stimulus onset (default: TRUE).
#' @param return_plot Logical indicating whether to return the ggplot object (TRUE) or print it (FALSE) (default: FALSE).
#' @return A ggplot object showing the spike raster.
#' @export
plot_raster <- function(
    nrn,
    plot_title = "Spike raster",
    zero_as_onset = TRUE,
    return_plot = FALSE
  ) {
    
    # Grab raster and time bounds
    spike.raster <- as.data.frame(nrn$fetch_spike_raster_R())
    total_time <- max(spike.raster$time) - min(spike.raster$time)
    dx <- total_time * 0.0025
    
    # Make plot
    title_size <- 14 
    axis_size <- 12 
    legend_size <- 10
    plt <- ggplot2::ggplot(spike.raster) +
      ggplot2::geom_segment(
        ggplot2::aes(x = .data$time - dx, xend = .data$time + dx, y = .data$trial, yend = .data$trial),
        color = "black",
        linewidth = 0.5
      ) +
      ggplot2::labs(
        title = make.plot.title(nrn, plot_title),
        x = "Time (ms)",
        y = "Trial"
      ) + 
      ggplot2::theme_minimal()
    if (zero_as_onset && min(spike.raster$time) <= 0) plt <- plt + 
      ggplot2::geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "darkblue") 
    plt <- plt + 
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
        axis.title = ggplot2::element_text(size = axis_size),
        axis.text = ggplot2::element_text(size = axis_size),
        legend.title = ggplot2::element_text(size = legend_size),
        legend.text = ggplot2::element_text(size = legend_size),
        legend.position = "bottom"
      )
    
    if (return_plot) {
      return(plt)
    } else {
      print(plt)
      return(invisible(NULL))
    }

  }

# S3/S4 methods for neuron objects ######################################################################################
#
# neuron objects are Rcpp Modules reference objects: formally S4 (class "Rcpp_neuron"), not
# plain S3. Prepending a class via class(x) <- c("neuron", class(x)) silently strips their
# S4-ness (R's own class<- warns "result will no longer be an S4 object"), which breaks the $
# accessor Rcpp modules rely on for method dispatch -- so these methods target "Rcpp_neuron"
# directly rather than a wrapper class. print()'s S3 dispatch only fires on an explicit
# print(nrn) call; auto-printing a neuron at the console (typing nrn by itself) goes through
# the S4 generic show() instead, hence the separate setMethod() below.

#' Print a neuron object
#'
#' @param x A neuron object.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly.
#' @export
print.Rcpp_neuron <- function(x, ...) {

    id_data <- x$fetch_id_data()
    has_trial_data <- nrow(x$fetch_trial_data_R()) > 0
    has_autocorr <- length(x$fetch_autocorr_R()) > 0

    header <- paste0("<neuron> #", id_data$id_num, " (", id_data$recording_name, ")")
    if (id_data$sim) header <- paste0(header, " [simulated]")
    cat(header, "\n")
    cat("  type:", id_data$type, "\n")
    for (field in c("genotype", "sex", "hemi", "region", "age")) {
      val <- id_data[[field]]
      if (!(val %in% c("not_provided", "not provided"))) cat("  ", field, ": ", val, "\n", sep = "")
    }
    cat("  trial data loaded:", has_trial_data, "\n")
    if (has_trial_data) {
      lambda <- x$fetch_lambda()
      cat("    lambda:    ", signif(unname(lambda["lambda"]), 4), id_data$unit_data, "per", id_data$unit_time, "\n")
      cat("    lambda_bin:", signif(unname(lambda["lambda_bin"]), 4), id_data$unit_data, "per bin\n")
    }
    cat("  autocorrelation computed:", has_autocorr, "\n")

    invisible(x)

  }

#' Summarize a neuron object
#'
#' @param object A neuron object.
#' @param ... Unused; present for S3 consistency.
#' @return A \code{summary.neuron} object, printed by its own \code{print} method.
#' @export
summary.Rcpp_neuron <- function(object, ...) {

    id_data <- object$fetch_id_data()
    has_trial_data <- nrow(object$fetch_trial_data_R()) > 0
    has_autocorr <- length(object$fetch_autocorr_R()) > 0
    has_fit <- length(object$fetch_autocorr_edf_R()) > 0

    out <- list(
      id_data = id_data,
      has_trial_data = has_trial_data,
      has_autocorr = has_autocorr,
      has_fit = has_fit,
      lambda = if (has_trial_data) object$fetch_lambda() else NULL,
      edf_parameters = if (has_fit) object$fetch_EDF_parameters() else NULL
    )
    class(out) <- "summary.neuron"
    out

  }

#' Print a neuron summary
#'
#' @param x A \code{summary.neuron} object, from \code{summary.neuron}.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly.
#' @export
print.summary.neuron <- function(x, ...) {

    id_data <- x$id_data
    cat("Neuron #", id_data$id_num, " (", id_data$recording_name, ")\n", sep = "")
    cat("  type:", id_data$type, "\n")
    cat("  simulated:", id_data$sim, "\n")
    for (field in c("genotype", "sex", "hemi", "region", "age")) {
      val <- id_data[[field]]
      if (!(val %in% c("not_provided", "not provided"))) cat("  ", field, ": ", val, "\n", sep = "")
    }
    cat("  unit_time:", id_data$unit_time, "\n")
    cat("  t_per_bin:", id_data$t_per_bin, id_data$unit_time, "\n")
    cat("  sample_rate:", id_data$sample_rate, id_data$unit_sample_rate, "\n")
    cat("  trial data loaded:", x$has_trial_data, "\n")
    if (x$has_trial_data) {
      cat("    lambda:    ", signif(unname(x$lambda["lambda"]), 4), id_data$unit_data, "per", id_data$unit_time, "\n")
      cat("    lambda_bin:", signif(unname(x$lambda["lambda_bin"]), 4), id_data$unit_data, "per bin\n")
    }
    cat("  autocorrelation computed:", x$has_autocorr, "\n")
    cat("  exponential decay fit:", x$has_fit, "\n")
    if (x$has_fit) {
      edf <- x$edf_parameters
      cat("    A:        ", signif(unname(edf["A"]), 4), "\n")
      cat("    tau:      ", signif(unname(edf["tau"]), 4), id_data$unit_time, "\n")
      cat("    bias_term:", signif(unname(edf["bias_term"]), 4), "\n")
    }

    invisible(x)

  }

#' Plot a neuron object
#'
#' Default \code{plot} method for neuron objects. Dispatches to \code{\link{plot_raster}} or \code{\link{plot_autocorrelation}} depending on \code{type}; defaults to the raster view.
#'
#' @param x A neuron object.
#' @param type One of "raster" (default) or "autocorrelation", specifying which view of the neuron to plot.
#' @param ... Additional arguments passed on to \code{\link{plot_raster}} or \code{\link{plot_autocorrelation}}.
#' @return Whatever \code{\link{plot_raster}}/\code{\link{plot_autocorrelation}} returns (a ggplot object if \code{return_plot = TRUE}, otherwise \code{NULL}, invisibly).
#' @export
plot.Rcpp_neuron <- function(x, type = c("raster", "autocorrelation"), ...) {
    type <- match.arg(type)
    if (type == "raster") {
      plot_raster(x, ...)
    } else {
      plot_autocorrelation(x, ...)
    }
  }

#' Helper function for testing assumption of autocorrelation processing
#'
#' The approach we're using is valid if autocorr (red lines) is below mean firing rate (blue line). If so, then the equation we're using to convert spiking autocorr to guassian sigma should work. Thin wrapper around \code{plot.neuronsDG_autocorr_table}, which does the actual plotting.
#'
#' @param autocor_results Data frame of autocorrelation results from the \code{process.autocorr} function.
#' @return Nothing; this function produces a base-R plot for visual inspection.
#' @export
test.sigma.assumption <- function(
    autocor_results
  ) {
    if (!inherits(autocor_results, "neuronsDG_autocorr_table")) {
      class(autocor_results) <- c("neuronsDG_autocorr_table", class(autocor_results))
    }
    plot(autocor_results)
  }

#' Plot an autocorrelation summary table
#'
#' \code{plot} method for the data frame returned by \code{\link{process.autocorr}}. Checks the assumption behind the sigma-fitting equation used to convert spiking autocorrelation to a dichotomized-Gaussian sigma: the approach is valid if autocorrelation (red lines) stays below mean firing rate (blue line).
#'
#' @param x A \code{neuronsDG_autocorr_table} object, i.e. the output of \code{\link{process.autocorr}}.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly. Produces a base-R plot as a side effect.
#' @export
plot.neuronsDG_autocorr_table <- function(x, ...) {

    x <- x[order(x$lambda_bin),]
    these_cols <- c(
      "lambda_bin",
      "max_autocorr",
      "mean_autocorr",
      "min_autocorr"
    )
    plot(
      x[,2],
      ylim = c(0, max(x[,these_cols])*1.1),
      type = "l",
      ylab = "value",
      xlab = "neuron (sorted by lambda_bin)",
      col = "blue"
    )
    lines(x$max_autocorr, col = "red4")
    lines(x$mean_autocorr, col = "red3")
    lines(x$min_autocorr, col = "red")

    invisible(x)

  }

#' Function to process autocorrelation of neuron list
#' 
#' First estimates autocorrelation of each neuron in the list, then fits an exponential decay function to the estimated autocorrelation, and finally collects and returns summary statistics for each neuron.
#' 
#' @param neuron_list An R list of neuron objects.
#' @param bin_count_action Method for counting spikes in each bin when computing autocorrelation; one of "boolean", "mean", or "sum" (default: "sum").
#' @param max_lag Maximum lag (in units of the trial data) to compute autocorrelation; if 0, uses the number of time bins in the neuron's trial data (default: 0).
#' @param A0 Initial guess for amplitude parameter of exponential decay function (default: 0.001).
#' @param tau0 Initial guess for time constant parameter of exponential decay function (default: 1.0).
#' @param ctol Convergence tolerance for fitting exponential decay function (default: 1e-8).
#' @param max_evals Maximum number of evaluations for fitting exponential decay function (default: 500).
#' @param check_autofiring_ratio Logical indicating whether to check the assumption that autocorrelation values are below the mean firing rate with \code{test.sigma.assumption} function (default: FALSE).
#' @param print_plots Logical indicating whether to print autocorrelation plots for each neuron (default: FALSE).
#' @param plot_time_cutoff Maximum lag (in bins) to display on the x-axis of plots (default: Inf).
#' @param use_raw Logical indicating whether to use raw autocorrelation (true) or standard centered and normalized correlation (false) (default: TRUE).
#' @return A data frame with one row per neuron and columns for lambda (mean spike rate per ms and per bin), amplitude (A), time constant (tau), bias term, first autocorrelation value, maximum autocorrelation value, mean autocorrelation value, and minimum autocorrelation value.
#' @export
process.autocorr <- function(
    neuron_list,
    bin_count_action = "sum", 
    max_lag = 0,
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500,
    check_autofiring_ratio = FALSE,
    print_plots = FALSE,
    plot_time_cutoff = Inf,
    use_raw = TRUE
  ) {
    
    # Array to hold results
    autocor_results <- c()
    
    # Find max lag if needed
    if (max_lag == 0) {
      trial_lengths <- rep(Inf, length(neuron_list))
      for (i in seq_along(neuron_list)) {
        nrn <- neuron_list[[i]]
        trial_lengths[i] <- nrow(nrn$fetch_trial_data_R())
      }
      max_lag <- min(trial_lengths)
    }
    
    # Loop through neurons in the list
    for (i in seq_along(neuron_list)) {
      
      # Grab neuron
      nrn <- neuron_list[[i]]
      
      # Set exponential decay function (EDF) parameters 
      nrn$set_edf_initials(A0, tau0)
      nrn$set_edf_termination(ctol, max_evals)
      
      # Compute autocorrelation
      nrn$compute_autocorrelation(bin_count_action, max_lag, use_raw)
      
      # Fit autocorrelation 
      nrn$fit_autocorrelation()
      
      # Fetch results 
      edf_results <- nrn$fetch_EDF_parameters()
      # ... get amplitude
      A <- edf_results["A"]
      # ... get time constant
      tau <- edf_results["tau"]
      # ... get bias term
      bias_term <- edf_results["bias_term"]
      # ... get lambda (mean spike rate, per ms and per bin)
      lambda <- nrn$fetch_lambda()
      # ... get series of EDF-fitted autocorrelation values
      autocorr_edf <- nrn$fetch_autocorr_edf_R()
      
      # Collect results and separate out the first predicted autocorrelation value
      autocor_results_r <- c(
        lambda, # vector of length two
        A,
        tau,
        bias_term,
        autocorr_edf[1],
        max(autocorr_edf[2:length(autocorr_edf)]),
        mean(autocorr_edf[2:length(autocorr_edf)]), 
        min(autocorr_edf[2:length(autocorr_edf)])
      )
      
      # Save results with those from other neurons
      autocor_results <- rbind(autocor_results, autocor_results_r)
      
      # Print plots 
      if (print_plots) {
        
        # Grab neuron ID info
        id_data <- nrn$fetch_id_data() 
        recording_name <- id_data$recording_name
        sim <- id_data$sim
        if (sim) recording_name <- paste(recording_name, "(sim)")
        
        print(plot_autocorrelation(nrn, "Est. autocorr", bias_term, plot_time_cutoff))
        
      }
      
    }
    
    autocor_results <- data.frame(names(neuron_list), autocor_results)
    
    rownames(autocor_results) <- NULL
    colnames(autocor_results) <- c(
      "cell",
      "lambda_ms",
      "lambda_bin",
      "A", 
      "tau", 
      "bias_term",
      "autocorr1",
      "max_autocorr",
      "mean_autocorr",
      "min_autocorr"
    )
    class(autocor_results) <- c("neuronsDG_autocorr_table", class(autocor_results))

    # Convert to data frame
    if (check_autofiring_ratio) {
      test.sigma.assumption(autocor_results)
    }

    return(autocor_results)

  }

#' Function to estimate autocorrelation parameters using dichotomized Gaussian simulations
#' 
#' This function performs the same estimate-and-fit procedure as \code{process.autocorr}, but does so multiple times for each neuron using simulated spike rasters generated from a dichotomized Gaussian model of that neuron. 
#' 
#' @param neuron_list An R list of neuron objects.
#' @param n_trials_per_sim Number of trials to simulate for each simulation (default: 300).
#' @param n_sims_per_neurons Number of simulations to run for each neuron (default: 100).
#' @param bin_count_action Method for counting spikes in each bin when computing autocorrelation; one of "boolean", "mean", or "sum" (default: "sum").
#' @param max_lag Maximum lag (in raw, unbinned units of each neuron's own trial data) to compute autocorrelation for both the original and each dichotomized-Gaussian-simulated spike raster; if 0, uses the shortest raw trial length across \code{neuron_list} (default: 0).
#' @param A0 Initial guess for amplitude parameter of exponential decay function (default: 0.001).
#' @param tau0 Initial guess for time constant parameter of exponential decay function (default: 1.0).
#' @param ctol Convergence tolerance for fitting exponential decay function (default: 1e-8).
#' @param max_evals Maximum number of evaluations for fitting exponential decay function (default: 500).
#' @param use_raw Logical indicating whether to use raw autocorrelation (true) or standard centered and normalized correlation (false) (default: TRUE).
#' @param seed Integer seed for the dichotomized-Gaussian simulation draws (default: 12345). These draws are made by a self-contained C++ random number generator, not R's own, so \code{set.seed()} has no effect on them -- reproducibility is controlled entirely by this argument. Each neuron in \code{neuron_list} gets its own seed, offset from this one so that neurons don't share a draw sequence.
#' @return A list containing a data frame of autocorrelation parameter estimates (one row per simulation), a data frame of neuron identifiers (one row per neuron), and the number of simulations run per neuron.
#' @export
estimate.autocorr.params <- function(
    neuron_list,
    n_trials_per_sim = 300,
    n_sims_per_neurons = 100,
    bin_count_action = "sum", # "boolean", "mean", or "sum"
    max_lag = 0,
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500,
    use_raw = TRUE,
    seed = 12345
  ) {
    
    # Find max lag if needed
    # ... max_lag must be in the same (raw, unbinned) time units as each neuron's
    #      trial data, matching what compute_autocorrelation() expects for the
    #      original (non-simulated) neuron; see neuron::estimate_autocorr_params,
    #      which separately converts to bin units for the simulated neurons.
    if (max_lag == 0) {
      trial_lengths <- rep(Inf, length(neuron_list))
      for (i in seq_along(neuron_list)) {
        nrn <- neuron_list[[i]]
        trial_lengths[i] <- nrow(nrn$fetch_trial_data_R())
      }
      max_lag <- min(trial_lengths)
    }

    # Pre-allocate matrix
    n_neurons <- length(neuron_list)
    dg_estimates <- matrix(NA_real_, nrow = n_sims_per_neurons * n_neurons, ncol = 9)
    
    # Initialize neuron id df
    neuron_id <- data.frame()
    
    # Loop through neurons in the list
    for (i in seq_along(neuron_list)) {
      
      # Get neuron 
      nrn <- neuron_list[[i]]
      
      # Get neuron id 
      nrn_id <- as.data.frame(lapply(nrn$fetch_id_data(), rep, times = 1))
      neuron_id <- rbind(neuron_id, nrn_id)
      
      # Run simulations to estimate value for this neuron and load into matrix
      row_idx_initial <- (i - 1) * n_sims_per_neurons + 1
      row_idx_final <- i * n_sims_per_neurons
      row_idx <- row_idx_initial:row_idx_final
      dg_estimates[row_idx,] <- nrn$estimate_autocorr_params(
        n_trials_per_sim,
        n_sims_per_neurons,
        max_lag,
        bin_count_action,
        A0,
        tau0,
        ctol,
        max_evals,
        use_raw,
        FALSE,
        as.integer(seed) + i # distinct, deterministic seed per neuron
      )
      
    }
    
    dg_estimates <- as.data.frame(dg_estimates)
    colnames(dg_estimates) <- c(
      "lambda_ms",
      "lambda_bin",
      "A", 
      "tau", 
      "bias_term",
      "autocorr1",
      "max_autocorr",
      "mean_autocorr",
      "min_autocorr"
    )
    
    result <- list(
      estimates = dg_estimates,
      neuron_id = neuron_id,
      n_sims_per_neurons = n_sims_per_neurons
    )
    class(result) <- c("neuronsDG_estimates", class(result))
    return(result)

  }

#' Print an autocorrelation-parameter estimate set
#'
#' @param x A \code{neuronsDG_estimates} object, from \code{\link{estimate.autocorr.params}}.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly.
#' @export
print.neuronsDG_estimates <- function(x, ...) {

    cat("<neuronsDG_estimates>\n")
    cat("  neurons:          ", nrow(x$neuron_id), "\n")
    cat("  sims per neuron:  ", x$n_sims_per_neurons, "\n")
    cat("  total simulations:", nrow(x$estimates), "\n")
    cat("  mean tau:         ", signif(mean(x$estimates$tau), 4), "\n")
    cat("  mean A:           ", signif(mean(x$estimates$A), 4), "\n")

    invisible(x)

  }

#' Function to analyze autocorrelation estimates with bootstrapping
#'
#' This function takes the output from \code{estimate.autocorr.params} and performs a bootstrap analysis to estimate the distribution of the mean time constant (tau) for different levels of specified covariates.
#'
#' @param ests A \code{neuronsDG_estimates} object (the output of \code{estimate.autocorr.params}), or a list of such objects.
#' @param covariate Character string or vector of character strings specifying the covariate(s) to analyze (e.g., "hemi").
#' @param n_bs Number of bootstrap resamples to perform (default: 1e4).
#' @return A \code{neuronsDG_analysis} object: a list with two elements:
#'  \describe{
#'    \item{resamples}{data frame with bootstrap resamples of the mean tau for each combination of covariate levels}
#'    \item{distribution_plot}{ggplot object showing the distribution of mean tau for each combination of covariate levels}
#'  }
#' @export
analyze.autocorr <- function(
    ests,
    covariate,
    n_bs = 1e4
  ) {
    
    # Check input types
    error_message <- "Input is not a valid estimate.autocorr.params output or list of such output."
    if (inherits(ests, "neuronsDG_estimates")) {
      if (!(all(covariate %in% names(ests$neuron_id)))) {
        stop("Covariate not found in neuron_id.")
      }
      # ... this is single input, wrap in list
      ests <- list(ests)
    } else if (is.list(ests) && !is.data.frame(ests)) {
      # ... check if list of proper inputs
      for (i in seq_along(ests)) {
        if (!inherits(ests[[i]], "neuronsDG_estimates")) {
          stop(error_message)
        }
        if (!(all(covariate %in% names(ests[[i]]$neuron_id)))) {
          stop("Covariate not found in neuron_id.")
        }
      }
    } else {
      stop(error_message)
    }
    
    # Combine batches into big data frame
    ests_all <- data.frame()
    for (i in seq_along(ests)) {
      
      # Grab estimate list 
      est <- ests[[i]]
      
      # Expand neuron_id 
      est$neuron_id <- est$neuron_id[rep(1:nrow(est$neuron_id), each = est$n_sims_per_neurons), ]
      
      # Combine with estimates 
      estimates_id <- cbind(est$estimates, est$neuron_id)
      
      # Add to big data frame
      ests_all <- rbind(ests_all, estimates_id)
      
    }
    
    # Get covariate levels and their interactions
    est_covariates <- list() 
    est_covariate_levels <- list()
    for (c in c(covariate)) {
      est_covariates[[c]] <- ests_all[,c]
      est_covariate_levels[[c]] <- unique(ests_all[,c])
    }
    cov_X <- do.call(expand.grid, list(est_covariate_levels, stringsAsFactors = FALSE))
    n_interX <- nrow(cov_X)
    
    # Run bootstraps 
    resamples <- matrix(NA, nrow = n_bs, ncol = n_interX)
    n_est <- nrow(ests_all)
    interX_names <- c()
    for (i in 1:n_interX) {
      
      lvl <- cov_X[i,]
      interX_names <- c(interX_names, paste0(lvl, collapse = "_"))
      lvl_mask <- rep(TRUE, n_est)
      
      if (length(lvl) == 1) {
        lvl_mask <- lvl_mask & ests_all[,c] == lvl
      } else {
        for (c in c(covariate)) {
          lvl_mask <- lvl_mask & ests_all[,c] == cov_X[i,c]
        }
      }
      
      if (any(lvl_mask)) {
        n_cells <- length(unique(ests_all$id_num))
        for (j in 1:n_bs) {
          # Net effect is to resample from each neuron's estimates in a way that 
          #   accounts for the uncertainty estimated by the DG. E.g., if there are n neurons, 
          #   each resample has a 1/n chance of drawing from a given neuron N; but, as N
          #   is represented by n_sims values from the DG simulations, the probability of a given value
          #   being drawn to represent N is determined by the simulations. 
          these_samples <- sample(ests_all[lvl_mask, "tau"], n_cells, replace = TRUE)
          resamples[j,i] <- mean(these_samples)
        }
      }
      
    }
    resamples <- as.data.frame(resamples)
    colnames(resamples) <- interX_names
    
    resamples <- resamples[, colSums(!is.na(resamples)) > 0]
    
    # Plot
    df <- data.frame(
      value = unlist(resamples, use.names = FALSE),
      group = rep(names(resamples), each = nrow(resamples))
    )
    
    # Plot
    title_size <- 14 
    axis_size <- 12 
    legend_size <- 10
    distribution_plot <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) +
      ggplot2::geom_density(alpha = 0.6) +
      #scale_y_continuous(transform = "log1p") +
      ggplot2::labs(title = "Expected Time Constant", x = "ms", y = "Density") +
      ggplot2::theme_minimal() + 
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
        plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
        axis.title = ggplot2::element_text(size = axis_size),
        axis.text = ggplot2::element_text(size = axis_size),
        legend.title = ggplot2::element_text(size = legend_size),
        legend.text = ggplot2::element_text(size = legend_size),
        legend.position = "bottom"
      ) 
    
    result <- list(
      resamples = resamples,
      distribution_plot = distribution_plot
    )
    class(result) <- c("neuronsDG_analysis", class(result))
    return(result)

  }

#' Print an autocorrelation bootstrap analysis
#'
#' @param x A \code{neuronsDG_analysis} object, from \code{\link{analyze.autocorr}}.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly.
#' @export
print.neuronsDG_analysis <- function(x, ...) {

    cat("<neuronsDG_analysis>\n")
    cat("  bootstrap resamples:", nrow(x$resamples), "\n")
    cat("  covariate groups:   ", paste(colnames(x$resamples), collapse = ", "), "\n")
    for (g in colnames(x$resamples)) {
      cat(sprintf("    %-20s mean tau = %s\n", g, signif(mean(x$resamples[[g]], na.rm = TRUE), 4)))
    }

    invisible(x)

  }

#' Plot an autocorrelation bootstrap analysis
#'
#' \code{plot} method for \code{neuronsDG_analysis} objects; displays the stored distribution plot.
#'
#' @param x A \code{neuronsDG_analysis} object, from \code{\link{analyze.autocorr}}.
#' @param ... Unused; present for S3 consistency.
#' @return \code{x}, invisibly. Prints \code{x$distribution_plot} as a side effect.
#' @export
plot.neuronsDG_analysis <- function(x, ...) {
    print(x$distribution_plot)
    invisible(x)
  }

#' Compute autocorrelation of single neuron object
#' 
#' This function computes the autocorrelation of a single neuron object using specified parameters. It's a wrapper around the \code{compute_autocorrelation} method of the neuron class.
#' 
#' @param nrn Neuron object for which to compute autocorrelation.
#' @param bin_count_action Method for counting spikes in each bin when computing autocorrelation; one of "boolean", "mean", or "sum" (default: "sum").
#' @param max_lag Maximum lag (in units of the trial data) to compute autocorrelation; if 0, uses the number of time bins in the neuron's trial data (default: 0).
#' @param use_raw Logical indicating whether to use raw autocorrelation (true) or standard centered and normalized correlation (false) (default: TRUE).
#' @return None. Modifies neuron in place.
#' @export
compute.autocorr <- function(
    nrn,
    bin_count_action = "sum",
    max_lag = 0,
    use_raw = TRUE
  ) {
    if (max_lag == 0) {
      max_lag <- nrow(nrn$fetch_trial_data_R())
    }
    nrn <- nrn$compute_autocorrelation(bin_count_action, max_lag, use_raw)
  }

#' Fit exponential decay function to autocorrelation of single neuron object
#' 
#' This function fits an exponential decay function to the autocorrelation of a single neuron object using specified parameters. It's a wrapper around the \code{fit_autocorrelation} method of the neuron class.
#'
#' @param nrn Neuron object for which to fit autocorrelation.
#' @param A0 Initial guess for amplitude parameter of exponential decay function (default: 0.001).
#' @param tau0 Initial guess for time constant parameter of exponential decay function (default: 1.0).
#' @param ctol Convergence tolerance for fitting exponential decay function (default: 1e-8).
#' @param max_evals Maximum number of evaluations for fitting exponential decay function (default: 500).
#' @return None. Modifies neuron in place.
#' @export
fit.edf.autocorr <- function(
    nrn,
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500
  ) {
    nrn$set_edf_initials(A0, tau0)
    nrn$set_edf_termination(ctol, max_evals)
    nrn <- nrn$fit_autocorrelation()
  }

