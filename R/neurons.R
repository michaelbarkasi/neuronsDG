
# Title Page ###########################################################################################################
# 
# Neurons: A framework for neuron modelling, simulation, and analysis
# 
# By Mike Barkasi
# GNU GPLv3: https://www.gnu.org/licenses/gpl-3.0.en.html
#   Copyright (c) 2025
#   
# Front matter #########################################################################################################

# Prelim code:
#   (none)

# Prelim comments and warnings: 
#   1. This script defines the R6 object class Neuron and related helper R6 classes, along with methods for them. 
#         Neuron objects and their methods are intended to be a flexible, common framework for integrating, 
#         analyzing, and simulating different kinds of neuron data. 

#library(Rcpp)
#library(RcppEigen)

# Source Cpp ###########################################################################################################

#sourceCpp("src/neuron.cpp")

.onLoad <- function(libname, pkgname) {
  Rcpp::loadModule("neuron", TRUE)
}

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

# Load neuron data from raster file
load.rasters.as.neurons <- function(
    raster_df,
    bin_size = 20, 
    sample_rt = 1e4,
    time_cutoff = Inf
  ) {
    
    # Input: 
    #   raster_df: Data frame (or file name to csv importable as such), each row a spike; 
    #                 must have columns: cell, time_in_ms, trial
    #                 optional columns: recording_name, hemisphere, genotype, sex, region, age
    # Output: 
    #   neuron_list: List of neuron objects, one per cell
    
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
      
      # Extract recording name and hemisphere
      recording_name <- "not_provided"
      if ("recording_name" %in% colnames(raster_df)) {
        recording_name <- unique(raster_df$recording_name[c_mask])
        if (length(recording_name) != 1) recording_name <- "not_provided"
      } 
      hemi <- "not_provided"
      if ("hemisphere" %in% colnames(raster_df)) {
        hemi <- unique(raster_df$hemisphere[c_mask])
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
      ", neuron ", recording_name, ", #", id_num, 
      ", ", genotype, ", ", sex, ", ", hemi, " hemi, ", region, ", ", age
    )
    
    return(plot_title)
    
  }

# Function to plot autocorrelation of neuron
plot.autocorrelation <- function(
    nrn,                      # neuron object
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
      bin = seq(1, length(autocorr) - 1, 1)
    )
    df_temp <- df_temp[df_temp$bin <= plot_time_cutoff,]
    
    # Make plot
    plt <- ggplot(df_temp) +
      geom_line(aes(x = bin, y = autocorrelation), color = "blue") +  
      geom_line(aes(x = bin, y = autocorrelation_fitted), color = "red") +  
      geom_hline(yintercept = bias_term, linewidth = 2, linetype = "dotted", color = "darkgray") +
      labs(
        title = plot_title,
        x = "Lag (bins)",
        y = "Autocorrelation"
      ) + 
      theme_minimal()
    plot(plt)
    
  }

# Function to plot spike raster of neuron
plot.raster <- function(
    neuron,                   # neuron object
    plot_title = "Spike raster",
    zero_as_onset = TRUE
  ) {
    
    # Grab raster and time bounds
    spike.raster <- as.data.frame(neuron$fetch_spike_raster_R())
    total_time <- max(spike.raster$time) - min(spike.raster$time)
    dx <- total_time * 0.0025
    
    # Make plot
    plt <- ggplot(spike.raster) +
      geom_segment(
        aes(x = time - dx, xend = time + dx, y = trial, yend = trial), 
        color = "black",
        linewidth = 0.5
      ) +
      labs(
        title = make.plot.title(neuron, plot_title),
        x = "Time (ms)",
        y = "Trial"
      ) + 
      theme_minimal()
    if (zero_as_onset && min(spike.raster$time) <= 0) plt <- plt + 
      geom_vline(xintercept = 0, linetype = "solid", linewidth = 1, color = "darkblue") 
    
    plot(plt)
    
  }

# Helper function for testing assumption of autocorrelation processing
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
    # Approach we're using is valid if autocorr (red lines) is below mean 
    # firing rate (blue line). If so, then the equation we're using to 
    # convert spiking autocorr to guassian sigma should work. 
    
  }

# Function to process autocorrelation of neuron list
process.autocorr <- function(
    neuron_list,              # An R list of neuron objects
    bin_count_action = "sum", # "boolean", "mean", or "sum" 
    A0 = 0.001,
    tau0 = 1.0,
    ctol = 1e-8,
    max_evals = 500,
    check_autofiring_ratio = FALSE,
    return_summary = FALSE,
    print_plots = FALSE,
    plot_time_cutoff = Inf
  ) {
    
    # Output: Dataframe with one row per neuron
    
    # Array to hold results
    autocor_results <- c()
    
    # Loop through neurons in the list
    for (nrn in neuron_list) {
      
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
        
        plot.autocorrelation(nrn, "Est. autocorr", bias_term, plot_time_cutoff)
        
      }
      
    }
    
    colnames(autocor_results) <- c(
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
      test.sigma.assumption(as.data.frame(autocor_results))
    }
    
    if (return_summary) return(as.data.frame(autocor_results))
    
  }

estimate.autocorr.params <- function(
    neuron_list,              # An R list of neuron objects
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
    
    # Loop through neurons in the list
    for (i in seq_along(neuron_list)) {
      
      # Get neuron 
      nrn <- neuron_list[[i]]
      
      # Get neuron id 
      nrn_id <- as.data.frame(lapply(nrn$fetch_id_data(), rep, times = 1))
      neuron_id <- rbind(neuron_id, nrn_id)
      
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
      
      # Load estimates into matrix
      row_idx_initial <- (i - 1) * n_sims_per_neurons + 1
      row_idx_final <- i * n_sims_per_neurons
      row_idx <- row_idx_initial:row_idx_final
      dg_estimates[row_idx,] <- dg_estimates_nrn
      
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