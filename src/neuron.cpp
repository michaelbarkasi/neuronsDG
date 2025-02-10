
// neuron.cpp
#include "neuron.h"

/*
 * Helper functions
 */

CharacterVector enum_prefix(std::string prefix, int n) {
  CharacterVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = prefix + std::to_string(i + 1);
  }
  return result;
}

/*
 * Main neuron class
 */

// Constructor
neuron::neuron(
  const int id_num, 
  const std::string recording_name, 
  const std::string type, 
  const std::string hemi,
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
    hemi(hemi), 
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
  if (unit_data == "spike") {
    infer_raster();
  }
} 

// Loading trial data (Rcpp matrix)
void neuron::load_trial_data_R(const NumericMatrix& td) {
  trial_data = as<MatrixXd>(td);
  if (unit_data == "spike") {
    infer_raster();
  }
} 

// Loading spike raster (Eigen matrix)
void neuron::load_spike_raster(const MatrixXd& sr) {spike_raster = sr;} 

// Loading spike raster (Rcpp matrix)
void neuron::load_spike_raster_R(const NumericMatrix& sr) {
  
  // Check input form
  if (unit_data != "spike") {Rcpp::stop("Spike raster can only be loaded if unit_data is 'spike'");}
  if (sr.ncol() != 2) {Rcpp::stop("Spike raster must have two columns");}
  CharacterVector sr_colnames = colnames(sr); 
  if (sr_colnames[0] != "time" || sr_colnames[1] != "trial") {Rcpp::stop("Spike raster must have columns named 'time' and 'trial'");}
  
  // Save spike raster
  spike_raster = as<MatrixXd>(sr);
  
  // Infer trial data
  infer_trial();
  
} 

// Fetching trial data (Eigen matrix)
MatrixXd neuron::fetch_trial_data() const {return trial_data;} 

// Fetching trial data (Rcpp matrix)
NumericMatrix neuron::fetch_trial_data_R() const {
  NumericMatrix trial_data_R = wrap(trial_data);
  colnames(trial_data_R) = enum_prefix("trial", trial_data_R.ncol());
  rownames(trial_data_R) = enum_prefix("time", trial_data_R.nrow());
  return trial_data_R;
} 

// Fetching spike raster (Eigen matrix)
MatrixXd neuron::fetch_spike_raster() const {return spike_raster;} 

// Fetching spike raster (Eigen matrix)
NumericMatrix neuron::fetch_spike_raster_R() const {
  NumericMatrix spike_raster_R = wrap(spike_raster);
  colnames(spike_raster_R) = CharacterVector::create("time", "trial");
  return spike_raster_R;
} 

// Infer trial data from spike raster
void neuron::infer_trial() {
  
  // Find needed dimensions
  int num_trials = static_cast<int>(spike_raster.col(1).maxCoeff());
  int max_time = static_cast<int>(spike_raster.col(0).maxCoeff());
  max_time = ((max_time + 9) / 10) * 10; // round up to nearest 10
  
  // Resize trial_data and set to zero
  trial_data.resize(max_time, num_trials);
  trial_data.setZero();
  
  // Fill trial data with recorded spikes 
  for (int r = 0; r < spike_raster.rows(); r++) {
    trial_data(static_cast<int>(spike_raster(r, 0)) - 1, static_cast<int>(spike_raster(r, 1)) - 1) += 1.0;
  }
  
}

// Infer spike raster from trial data
void neuron::infer_raster() {
  
  int num_spikes = static_cast<int>(trial_data.sum());
  spike_raster.resize(num_spikes, 2);
  int ctr = 0;
  
  for (int r = 0; r < trial_data.rows(); r++) {
    for (int c = 0; c < trial_data.cols(); c++) {
      if (trial_data(r, c) > 0.0) {
        spike_raster(ctr, 0) = static_cast<double>(r + 1);
        spike_raster(ctr, 1) = static_cast<double>(c + 1);
        ctr++;
      }
    }
  }
  
}

List neuron::fetch_id_data() const {
  return List::create(
    _["id_num"] = id_num,
    _["recording_name"] = recording_name,
    _["type"] = type,
    _["hemi"] = hemi,
    _["sim"] = sim,
    _["unit_time"] = unit_time,
    _["unit_sample_rate"] = unit_sample_rate,
    _["unit_data"] = unit_data,
    _["t_per_bin"] = t_per_bin,
    _["sample_rate"] = sample_rate
  );
  
}

/*
 * Member function implementations, fetch analysis results
 */

// Fetching autocorrelation (Rcpp vector)
VectorXd neuron::fetch_autocorr() const {return autocorr;}

// Fetching autocorrelation (Rcpp vector)
NumericVector neuron::fetch_autocorr_R() const {return wrap(autocorr);}

/*
 * Member function implementations, analysis
 */

// Compute autocorrelation of trial data
void neuron::compute_autocorrelation(
    const std::string& bin_count_action
) {
  
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
  
  // Check binning action
  if (bin_count_action != "sum" && bin_count_action != "boolean" && bin_count_action != "mean") {
    Rcpp::stop("bin_count_action must be 'sum', 'boolean', or 'mean'");
  }
  
  // Initialize matrix for data manipulation 
  const int T_n = max_lag;
  const int N_trial = trial_data.cols();
  MatrixXd data(T_n, N_trial);
  
  // Collapse bins if needed
  if (convert_to_bins) {
    
    for (int n = 0; n < N_trial; n++) {
      for (int b = 0; b < T_n; b++) {
        double bin_count = trial_data(seq(b*(int)t_per_bin, (b*(int)t_per_bin) + t_per_bin - 1), n).sum();
        if (bin_count_action == "mean") {
          data(b, n) = bin_count/(double)t_per_bin;
        } else if (bin_count_action == "sum") {
          data(b, n) = bin_count;
        } else if (bin_count_action == "boolean") {
          data(b, n) = (bin_count > 0.0) ? 1.0 : 0.0;
        }
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
  .constructor<int, std::string, std::string, std::string, bool, std::string, std::string, std::string, double, double>()
  .method("load_trial_data_R", &neuron::load_trial_data_R)
  .method("load_spike_raster_R", &neuron::load_spike_raster_R)
  .method("fetch_trial_data_R", &neuron::fetch_trial_data_R)
  .method("fetch_spike_raster_R", &neuron::fetch_spike_raster_R)
  .method("fetch_autocorr_R", &neuron::fetch_autocorr_R)
  .method("fetch_id_data", &neuron::fetch_id_data)
  .method("compute_autocorrelation", &neuron::compute_autocorrelation);
}

