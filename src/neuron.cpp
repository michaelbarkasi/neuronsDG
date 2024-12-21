
// neuron.cpp
#include "neuron.h"

// Constructor
neuron::neuron(
  const int id_num, 
  const std::string recording_name, 
  const std::string type, 
  const bool sim, 
  const std::string unit_time, 
  const std::string unit_sample_rate, 
  const std::string unit_data, 
  const double t_per_bin, 
  const double sample_rate
) 
  : id_num(id_num), 
    recording_name(recording_name), 
    type(type), 
    sim(sim), 
    unit_time(unit_time), 
    unit_sample_rate(unit_sample_rate), 
    unit_data(unit_data), 
    t_per_bin(t_per_bin), 
    sample_rate(sample_rate)
{ 
  // No initialization operations
}

/*
 * Member function implementations, basic data handling
 */

// Loading trial data (Eigen matrix)
void neuron::load_trial_data(const MatrixXd& td) {
  trial_data = td;
} 

// Loading trial data (Rcpp matrix)
void neuron::load_trial_data_R(const NumericMatrix& td) {
  trial_data = as<MatrixXd>(td);
} 

// Loading spike raster (Eigen matrix)
void neuron::load_spike_raster(const MatrixXd& sr) {
  spike_raster = sr;
} 

// Loading spike raster (Rcpp matrix)
void neuron::load_spike_raster_R(const NumericMatrix& sr) {
  spike_raster = as<MatrixXd>(sr);
} 

// Fetching trial data (Eigen matrix)
MatrixXd neuron::fetch_trial_data() const {
  return trial_data;
} 

// Fetching trial data (Rcpp matrix)
NumericMatrix neuron::fetch_trial_data_R() const {
  return wrap(trial_data);
} 

// Fetching spike raster (Eigen matrix)
MatrixXd neuron::fetch_spike_raster() const {
  return spike_raster;
} 

// Fetching spike raster (Eigen matrix)
NumericMatrix neuron::fetch_spike_raster_R() const {
  return wrap(spike_raster);
} 

/*
 * Member function implementations, fetch analysis results
 */

// Fetching autocorrelation (Rcpp vector)
VectorXd neuron::fetch_autocorr() const {
  return autocorr;
}

// Fetching autocorrelation (Rcpp vector)
NumericVector neuron::fetch_autocorr_R() const {
  return wrap(autocorr);
}

/*
 * Member function implementations, analysis
 */

// Compute autocorrelation of trial data
void neuron::compute_autocorrelation() {
  
  /*
   * Based on Neophytou et al 2022, "Differences in temporal processing ... "
   *   https://doi.org/10.1371/journal.pbio.3001803
   *   equation: x_n(t) = trial_matrix[t,n]
   */
  
  // Check time unit of the trial matrix
  bool convert_to_bins = false;
  int max_lag = trial_data.rows();
  if (unit_time != "bin") {
    convert_to_bins = true;
    max_lag = (int) max_lag/t_per_bin;
  }
  
  // Initialize matrix for data manipulation 
  const int T_n = max_lag;
  const int N_trial = trial_data.cols();
  MatrixXd data(T_n, N_trial);
  
  // Collapse bins if needed
  if (convert_to_bins) {
    
    for (int n = 0; n < N_trial; n++) {
      for (int b = 0; b < T_n; b++) {
        data(b, n) = trial_data(seq(b*(int)t_per_bin, (b*(int)t_per_bin) + t_per_bin - 1), n).sum()/(double)t_per_bin;
      }
    }
    
  } else {
    data = trial_data;
  }
  
  // Pad the trial matrix with zeros
  MatrixXd data_padded(2 * data.rows(), data.cols());
  data_padded.topRows(data.rows()) = data;
  data_padded.bottomRows(data.rows()).setZero();
  
  // (Re)initialize autocorr
  autocorr.resize(max_lag);
  autocorr.setZero();
  
  // Compute autocorrelation for all possible lags
  for (int lag = 0; lag < max_lag; lag++) {
    
    double normalization_term = 1.0/((double)T_n - (double)lag + 1.0);
    // Sum over trials
    for (int n = 0; n < N_trial; n++) {
      autocorr(lag) += normalization_term * data_padded(seq(0 + lag, T_n + lag - 1),n).dot(data.col(n));
    }
    // Find mean proportion of co-active bins per trial
    autocorr(lag) *= (1.0/N_trial);
    
  }
  
}

/*
 * RCPP_MODULE to expose class to R and function to initialize neuron
 */

RCPP_MODULE(neuron) {
  class_<neuron>("neuron")
  .constructor<int, std::string, std::string, bool, std::string, std::string, std::string, double, double>()
  .method("load_trial_data_R", &neuron::load_trial_data_R)
  .method("load_spike_raster_R", &neuron::load_spike_raster_R)
  .method("fetch_trial_data_R", &neuron::fetch_trial_data_R)
  .method("fetch_spike_raster_R", &neuron::fetch_spike_raster_R)
  .method("fetch_autocorr_R", &neuron::fetch_autocorr_R)
  .method("compute_autocorrelation", &neuron::compute_autocorrelation);
}

