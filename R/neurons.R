
# Intro ################################################################################################################

# Neurons: A framework for neuron modelling, simulation, and analysis

# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2025

#' @useDynLib neurons, .registration = TRUE
#' @import Rcpp
#' @import RcppEigen 
#' @import ggplot2
#' @import parallel
NULL

.onLoad <- function(libname, pkgname) {
    Rcpp::loadModule("neuron", TRUE)
  }

core_num <- parallel::detectCores() - 2 # Set to 1 to disable parallel processing (e.g., on windows)

runparallel <- function(
    input, 
    funct, 
    ...
  ) {
    results <- parallel::mclapply( 
      X = input, FUN = funct, ...,
      mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = core_num,
      mc.cleanup = TRUE, mc.allow.recursive = TRUE 
      )
    return(results)
  }

# Define functions #####################################################################################################

#' Initialize neuron object
#' 
#' This function initializes a new neuron object with specified parameters. The neuron objects are the main tools for running autocorrelation analyses of spike data.
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
#' @param t_per_bin Time (in above units) per bin, e.g., 1 ms per bin (default: 1.0).
#' @param sample_rate Sample rate (in above units), e.g., 1e4 Hz (default: 1e4).
#' @return A new neuron object.
#' @export
new_neuron <- function(
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
    t_per_bin = 1.0, 
    sample_rate = 1e4
  ) {
    neuron <- new(
      neuron, 
      id_num, recording_name, type, genotype, sex, hemi, region, age, sim, unit_time, unit_sample_rate, unit_data, t_per_bin, sample_rate
    )
    return(neuron)
  }

#' Load neuron data from raster file
#' 
#' This function loads spike raster data from a data frame or CSV file and converts each unique cell into a neuron object. The raster data frame must contain columns for cell identifier, spike time in milliseconds, and trial number. Optional metadata columns can also be included.
#' 
#' @param raster_df Data frame (or file name to csv importable as such), each row a spike; must have columns: cell, time_in_ms, trial; optional columns: recording_name, hemisphere, genotype, sex, region, age.
#' @param bin_size Size of time bins in milliseconds (default: 20).
#' @param sample_rt Sample rate in the default unit for neuron objects, Hz (default: 1e4).
#' @param time_cutoff Maximum time (in ms) to include spikes; spikes occurring after this time will be excluded (default: Inf).
#' @return A list of neuron objects, one per unique cell in the raster data.
#' @export
load.rasters.as.neurons <- function(
    raster_df,
    bin_size = 20, 
    sample_rt = 1e4,
    time_cutoff = Inf
  ) {
    
    # Load data
    if (class(raster_df) == "character") {
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
      new_neuron <- new_neuron(
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
      new_neuron$load_spike_raster_R(raster)
      
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

#' Function to plot autocorrelation of neuron
#' 
#' Wrapper function to fetch and plot autocorrelation information from a neuron object. Assumes the autocorrelation has already been computed. 
#' 
#' @name plot.autocorrelation 
#' @rdname plot-autocorrelation
#' @usage plot.autocorrelation(nrn, plot_title = "Est. autocorr", bias_term = 0, plot_time_cutoff = Inf)
#' @param nrn Neuron object for which to plot autocorrelation.
#' @param plot_title Title for the plot (default: "Est. autocorr").
#' @param bias_term Bias term to plot as a horizontal line (default: 0).
#' @param plot_time_cutoff Maximum lag (in bins) to display on the x-axis (default: Inf).
#' @return A ggplot object showing the estimated and fitted autocorrelation.
#' @export
plot.autocorrelation <- function(
    nrn,
    plot_title = "Est. autocorr",
    bias_term = 0,
    plot_time_cutoff = Inf
  ) {
    
    # Make title
    plot_title <- make.plot.title(nrn, plot_title)
    
    # Fetch estimated and fitted autocorrelation 
    autocorr <- nrn$fetch_autocorr_R()
    autocorr_edf <- nrn$fetch_autocorr_edf_R()
   
    # Make data frame for plotting 
    df_temp <- data.frame(
      autocorrelation = autocorr[2:length(autocorr)],
      autocorrelation_fitted = autocorr_edf[2:length(autocorr_edf)],
      bin = seq(from = 1, to = length(autocorr) - 1, by = 1)
    )
    df_temp <- df_temp[df_temp$bin <= plot_time_cutoff,]
    
    # Make plot
    plt <- ggplot2::ggplot(df_temp) +
      ggplot2::geom_line(ggplot2::aes(x = bin, y = autocorrelation), color = "blue") +  
      ggplot2::geom_line(ggplot2::aes(x = bin, y = autocorrelation_fitted), color = "red") +  
      ggplot2::geom_hline(yintercept = bias_term, linewidth = 2, linetype = "dotted", color = "darkgray") +
      ggplot2::labs(
        title = plot_title,
        x = "Lag (bins)",
        y = "Autocorrelation"
      ) + 
      ggplot2::theme_minimal() + 
      ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
      )
    return(plt)
    
  }

#' Function to plot spike raster of neuron
#' 
#' Wrapper function to fetch and plot spike raster from a neuron object. Assumes the spike raster has already been loaded.
#' 
#' @name plot.raster
#' @rdname plot-raster
#' @usage plot.raster(nrn, plot_title = "Spike raster", zero_as_onset = TRUE)
#' @param nrn Neuron object for which to plot spike raster.
#' @param plot_title Title for the plot (default: "Spike raster").
#' @param zero_as_onset Logical indicating whether to plot a vertical line at time zero to indicate stimulus onset (default: TRUE).
#' @return A ggplot object showing the spike raster.
#' @export
plot.raster <- function(
    nrn,  
    plot_title = "Spike raster",
    zero_as_onset = TRUE
  ) {
    
    # Grab raster and time bounds
    spike.raster <- as.data.frame(nrn$fetch_spike_raster_R())
    total_time <- max(spike.raster$time) - min(spike.raster$time)
    dx <- total_time * 0.0025
    
    # Make plot
    plt <- ggplot2::ggplot(spike.raster) +
      ggplot2::geom_segment(
        ggplot2::aes(x = time - dx, xend = time + dx, y = trial, yend = trial), 
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
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
      )
    
    return(plt)
    
  }

#' Helper function for testing assumption of autocorrelation processing
#' 
#' The approach we're using is valid if autocorr (red lines) is below mean firing rate (blue line). If so, then the equation we're using to convert spiking autocorr to guassian sigma should work. 
#' 
#' @param autocor_results Data frame of autocorrelation results from the \code{process.autocorr} function.
#' @return Nothing; this function produces a base-R plot for visual inspection.
#' @export
test.sigma.assumption <- function(
    autocor_results
  ) {
    
    autocor_results <- autocor_results[order(autocor_results$lambda_bin),]
    these_cols <- c(
      "lambda_bin",
      "max_autocorr",
      "mean_autocorr",
      "min_autocorr"
    )
    plot(
      autocor_results[,2], 
      ylim = c(0, max(autocor_results[,these_cols])*1.1),
      type = "l",
      ylab = "value",
      xlab = "neuron (sorted by lambda_bin)",
      col = "blue"
    )
    lines(autocor_results$max_autocorr, col = "red4")
    lines(autocor_results$mean_autocorr, col = "red3")
    lines(autocor_results$min_autocorr, col = "red")
    
  }

#' Function to process autocorrelation of neuron list
#' 
#' First estimates autocorrelation of each neuron in the list, then fits an exponential decay function to the estimated autocorrelation, and finally collects and returns summary statistics for each neuron.
#' 
#' @param neuron_list An R list of neuron objects.
#' @param bin_count_action Method for counting spikes in each bin when computing autocorrelation; one of "boolean", "mean", or "sum" (default: "sum").
#' @param A0 Initial guess for amplitude parameter of exponential decay function (default: 0.001).
#' @param tau0 Initial guess for time constant parameter of exponential decay function (default: 1.0).
#' @param ctol Convergence tolerance for fitting exponential decay function (default: 1e-8).
#' @param max_evals Maximum number of evaluations for fitting exponential decay function (default: 500).
#' @param check_autofiring_ratio Logical indicating whether to check the assumption that autocorrelation values are below the mean firing rate with \code{test.sigma.assumption} function (default: FALSE).
#' @param print_plots Logical indicating whether to print autocorrelation plots for each neuron (default: FALSE).
#' @param plot_time_cutoff Maximum lag (in bins) to display on the x-axis of plots (default: Inf).
#' @return A data frame with one row per neuron and columns for lambda (mean spike rate per ms and per bin), amplitude (A), time constant (tau), bias term, first autocorrelation value, maximum autocorrelation value, mean autocorrelation value, and minimum autocorrelation value.
#' @export
process.autocorr <- function(
    neuron_list,
    bin_count_action = "sum", 
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500,
    check_autofiring_ratio = FALSE,
    print_plots = FALSE,
    plot_time_cutoff = Inf
  ) {
    
    # Array to hold results
    autocor_results <- c()
    
    # Loop through neurons in the list
    for (i in seq_along(neuron_list)) {
      
      nrn <- neuron_list[[i]]
      
      # Set exponential decay function (EDF) parameters 
      nrn$set_edf_initials(A0, tau0)
      nrn$set_edf_termination(ctol, max_evals)
      
      # Compute autocorrelation
      nrn$compute_autocorrelation(bin_count_action)
      
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
        
        plot(plot.autocorrelation(nrn, "Est. autocorr", bias_term, plot_time_cutoff))
        
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
#' @param A0 Initial guess for amplitude parameter of exponential decay function (default: 0.001).
#' @param tau0 Initial guess for time constant parameter of exponential decay function (default: 1.0).
#' @param ctol Convergence tolerance for fitting exponential decay function (default: 1e-8).
#' @param max_evals Maximum number of evaluations for fitting exponential decay function (default: 500).
#' @param verbose Logical indicating whether to print progress messages (default: FALSE).
#' @return A list containing a data frame of autocorrelation parameter estimates (one row per simulation), a data frame of neuron identifiers (one row per neuron), and the number of simulations run per neuron.
#' @export
estimate.autocorr.params <- function(
    neuron_list,
    n_trials_per_sim = 300, 
    n_sims_per_neurons = 100, 
    bin_count_action = "sum", # "boolean", "mean", or "sum" 
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500,
    verbose = FALSE
  ) {
    
    # Pre-allocate matrix 
    n_neurons <- length(neuron_list)
    dg_estimates <- matrix(NA_real_, nrow = n_sims_per_neurons * n_neurons, ncol = 9)
    
    # Initialize neuron id df
    neuron_id <- data.frame()
    
    est_autocorr_params <- function(i) {
      # Get neuron 
      nrn <- neuron_list[[i]]
      # Run simulations to estimate value for this neuron
      dg_estimates_nrn <- nrn$estimate_autocorr_params(
        n_trials_per_sim,
        n_sims_per_neurons,
        bin_count_action,
        A0,
        tau0,
        ctol,
        max_evals,
        verbose
      )
      return(dg_estimates_nrn)
    }
    
    dg_estimates_nrn_list <- runparallel(seq_along(neuron_list), est_autocorr_params)
    
    # Loop through neurons in the list
    for (i in seq_along(neuron_list)) {
      
      # Get neuron 
      nrn <- neuron_list[[i]]
      
      # Get neuron id 
      nrn_id <- as.data.frame(lapply(nrn$fetch_id_data(), rep, times = 1))
      neuron_id <- rbind(neuron_id, nrn_id)
      
      # Load estimates into matrix
      row_idx_initial <- (i - 1) * n_sims_per_neurons + 1
      row_idx_final <- i * n_sims_per_neurons
      row_idx <- row_idx_initial:row_idx_final
      dg_estimates[row_idx,] <- dg_estimates_nrn_list[[i]]
      
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
    
    return(list(
      estimates = dg_estimates,
      neuron_id = neuron_id,
      n_sims_per_neurons = n_sims_per_neurons
    ))
    
  }

#' Function to analyze autocorrelation estimates with bootstrapping
#' 
#' This function takes the output from \code{estimate.autocorr.params} and performs a bootstrap analysis to estimate the distribution of the mean time constant (tau) for different levels of specified covariates.
#' 
#' @param ests Output from \code{estimate.autocorr.params}, or a list of such outputs.
#' @param covariate Character string or vector of character strings specifying the covariate(s) to analyze (e.g., "hemi").
#' @param n_bs Number of bootstrap resamples to perform (default: 1e4).
#' @return A list with two elements: 
#'  \describe{
#'    \item{resamples}{data frame with bootstrap resamples of the mean tau for each combination of covariate levels}
#'    \item{distribution_plot}{ggplot object showing the distribution of mean tau for each combination of covariate levels}
#'  }
analyze.autocorr <- function(
    ests,
    covariate,
    n_bs = 1e4
  ) {
    
    # Check input types
    error_message <- "Input is not a valid estimate.autocorr.params output or list of such output."
    if (class(ests) == "list") {
      if ("n_sims_per_neurons" %in% names(ests)) {
        if (!(all(covariate %in% names(ests$neuron_id)))) {
          stop("Covariate not found in neuron_id.")
        }
        # ... this is single input, wrap in list
        ests <- list(ests)
      } else {
        # ... check if list of proper inputs 
        for (i in seq_along(ests)) {
          if (!("n_sims_per_neurons" %in% names(ests[[i]]))) {
            stop(error_message)
          }
          if (!(all(covariate %in% names(ests[[i]]$neuron_id)))) {
            stop("Covariate not found in neuron_id.")
          }
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
      est$estimates <- cbind(est$estimates, est$neuron_id)
      
      # Add to big data frame
      ests_all <- rbind(ests_all, est$estimates)
      
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
          lvl_mask <- lvl_mask & ests_all[,c] == lvl[[c]]
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
          resamples[j,i] <- mean(sample(ests_all[lvl_mask, "tau"], n_cells, replace = TRUE))
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
    title_size <- 20 
    axis_size <- 12 
    legend_size <- 10
    distribution_plot <- ggplot2::ggplot(df, ggplot2::aes(x = value, fill = group)) +
      ggplot2::geom_density(alpha = 0.6) +
      #scale_y_continuous(transform = "log1p") +
      ggplot2::labs(title = "Expected Time Constant", x = "ms", y = "Density") +
      ggplot2::theme_minimal() + 
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = title_size),
        axis.title = ggplot2::element_text(size = axis_size),
        axis.text = ggplot2::element_text(size = axis_size),
        legend.title = ggplot2::element_text(size = legend_size),
        legend.text = ggplot2::element_text(size = legend_size),
        legend.position = "bottom"
      ) + ggplot2::theme(
        panel.background = ggplot2::element_rect(fill = "white", colour = NA),
        plot.background  = ggplot2::element_rect(fill = "white", colour = NA)
      )
    
    return(
      list(
        resamples = resamples,
        distribution_plot = distribution_plot
      )
    )
    
  }
