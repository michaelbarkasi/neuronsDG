
// neuron.cpp
#include "neuron.h"

/*
 * ***********************************************************************************
 * Helper functions
 */

// Build sequence of numbered string prefixes
CharacterVector enum_prefix(std::string prefix, int n) {
  CharacterVector result(n);
  for (int i = 0; i < n; ++i) {
    result[i] = prefix + std::to_string(i + 1);
  }
  return result;
}

// Convert to std::vector with doubles 
std::vector<double> to_dVec(
    const VectorXd& vec
  ) {
    std::vector<double> dVec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      dVec[i] = vec(i);
    }
    return dVec;
  }
// ... overload
std::vector<double> to_dVec(
    const NumericVector& vec
  ) {
    return Rcpp::as<std::vector<double>>(vec);
  }

// Convert to Eigen vector with doubles
VectorXd to_eVec(
    const std::vector<double>& vec
  ) {
    VectorXd eVec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      eVec(i) = vec[i];
    }
    return eVec;
  }

// Convert to NumericVector 
NumericVector to_NumVec(
    const VectorXd& vec
  ) {
    NumericVector num_vec(vec.size());
    for (int i = 0; i < vec.size(); i++) {
      num_vec(i) = vec(i);
    }
    return num_vec;
  }

// Convert to Eigen matrix with doubles
MatrixXd to_eMat(
    const NumericMatrix& X
  ) {
    int Xnrow = X.nrow();
    int Xncol = X.ncol();
    MatrixXd M = MatrixXd(Xnrow, Xncol);
    for (int i = 0; i < Xnrow; i++) {
      for (int j = 0; j < Xncol; j++) {
        M(i, j) = X(i, j);
      }
    }
    return M;
  }

// Convert to NumericMatrix
NumericMatrix to_NumMat(
    const MatrixXd& M
  ) {
    int M_nrow = M.rows();
    int M_ncol = M.cols();
    NumericMatrix X(M_nrow, M_ncol);
    for (int i = 0; i < M_nrow; i++) {
      for (int j = 0; j < M_ncol; j++) {
        X(i, j) = M(i, j);
      }
    }
    return X;
  }

// EDF model 
double EDF_autocorr(
    const double& lag,
    const double& A, 
    const double& tau, 
    const double& bias_term,
    const int& return_grad
  ) {
    if (return_grad == 0) { // No gradient, return function output
      return A * exp(-lag/tau) + bias_term;
    } else if (return_grad == 1) { // gradient wrt A
      return exp(-lag/tau);
    } else if (return_grad == 2) { // gradient wrt tau
      return A * exp(-lag/tau) * (lag/(tau*tau));
    } else {
      return 0.0;
    }
  }

// Multivariate normal CDF
double mvnorm_cdf(
    const NumericVector& lower, 
    const NumericVector& upper, 
    const NumericVector& mean, 
    const NumericMatrix& sigma  // covariance matrix of dimension n less than 1000
  ) {
    
    if (sigma.nrow() != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (sigma.nrow() >= 1000) {Rcpp::stop("Covariance matrix must be less than 1000x1000");}
    if (lower.size() != upper.size()) {Rcpp::stop("Lower and upper bounds must be same length");}
    if (mean.size() != sigma.nrow() || mean.size() != lower.size()) {Rcpp::stop("Bounds and mean vectors must have the same length as sigma diagonal");}
    
    Function pmvnorm("pmvnorm"); // Calls mvtnorm::pmvnorm from R
    // ... uses Genz algorithm
    
    return as<double>(
      pmvnorm(
        Named("lower") = lower, 
        Named("upper") = upper, 
        Named("mean") = mean, 
        Named("sigma") = sigma
      )
    );
    
  }

// Normal CDF inverse
double norm_cdf(
    const double& x,
    const double& mu,
    const double& sd,
    const bool& inverse
  ) {
    using boost::math::normal; 
    normal standard_normal(mu, sd);
    if (inverse) {
      return boost::math::quantile(standard_normal, x);
    } else {
      return boost::math::cdf(standard_normal, x);
    }
  }

// Multivariate normal random number generator
NumericMatrix mvnorm_random(
    int n,                   // Number of points to generate
    NumericVector mu,        // Mean vector, length determines dimension
    NumericMatrix sigma      // Covariance matrix, square, same dimension as mu
  ) {
   
    if (sigma.nrow() != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (sigma.nrow() != mu.size()) {Rcpp::stop("Mean vector must have same length as sigma diagonal");}
    
    Function mvrnorm("mvrnorm", Environment::namespace_env("MASS"));
    
    // Generate random points
    NumericMatrix X = as<NumericMatrix>(
      mvrnorm(
        Named("n") = n, 
        Named("mu") = mu, 
        Named("Sigma") = sigma,
        Named("tol") = 1e-6 // default is 1e-6
      )
    );
    
    return X;
    
  } 

// Function to create Toeplitz matrix
NumericMatrix toeplitz(
    const std::vector<double>& first_col, 
    const std::vector<double>& first_row
  ) {
    
    int rows = first_col.size();
    int cols = first_row.size();
    NumericMatrix T(rows, cols);
    
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        T(i, j) = (i >= j) ? first_col[i - j] : first_row[j - i];
      }
    }
    
    return T;
    
  }

// For estimating sigma for dichotomized Gaussian simulation
NumericVector dichot_gauss_sigma_formula(
    const double& threshold,      // threshold for dichotomization
    const NumericVector& cov,     // desired covarance after dichotomization
    const NumericMatrix& sigma    // covariance matrix
  ) {
   
    // Check dimension
    int dim = sigma.nrow();
    if (dim != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (dim != cov.size()) {Rcpp::stop("Covariance vector must have the same length as sigma diagonal");}
    
    // Find probability of a point being below the threshold along both dimensions
    NumericVector lower = Rcpp::rep(-1e100, dim);
    NumericVector upper = Rcpp::rep(threshold, dim); 
    double Phi2 = mvnorm_cdf(
      lower,
      upper, 
      0.0, 
      sigma
    );
    
    // Find probability of a point being below the threshold along one dimension
    double Phi = norm_cdf(
      threshold, 
      0.0,   // mean
      1,     // sd
      false  // return cdf
    );
    
    // desired sigma will be the one which sends all elements to zero
    NumericVector output(dim);
    double Phi2PhiPhi = Phi2 * Phi * Phi;
    for (int i = 0; i < dim; i++) {
      output[i] = cov[i] - Phi2PhiPhi;
    }
    return output;
    
  }

/*
 * ***********************************************************************************
 * Main neuron class
 */

// Constructor
neuron::neuron(
  const int id_num, 
  const std::string recording_name, 
  const std::string type, 
  const std::string hemi,
  bool sim, 
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
 * ***********************************************************************************
 * Member function implementations, adjust settings
 */

void neuron::set_edf_initials(
    double a0, 
    double t0
  ) {
    A0 = a0;
    tau0 = t0;
  };

void neuron::set_edf_termination(
    double ct, 
    int me
  ) {
    ctol = ct;
    max_evals = me;
  };

/*
 * ***********************************************************************************
 * Member function implementations, basic data handling
 */

// Loading trial data (Eigen matrix)
void neuron::load_trial_data(const MatrixXd& td) {
    
    // Set trial data
    trial_data = td;
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    
    // Make spike raaster
    if (unit_data == "spike") {
      infer_raster();
    }
    
  } 

// Loading trial data (Rcpp matrix)
void neuron::load_trial_data_R(const NumericMatrix& td) {
    
    // Set trial data
    trial_data = as<MatrixXd>(td);
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    
    // Make spike raster
    if (unit_data == "spike") {
      infer_raster();
    }
    
  } 

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
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    
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

// Method to return fields with neuron ID information
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
 * ***********************************************************************************
 * Member function implementations, fetch analysis results
 */

// Fetching autocorrelation 
VectorXd neuron::fetch_autocorr() const {return autocorr;}

// Fetching autocorrelation 
NumericVector neuron::fetch_autocorr_R() const {return wrap(autocorr);}
NumericVector neuron::fetch_autocorr_edf_R() const {return wrap(autocorr_edf);}

// Fetch fitted EDF parameters 
NumericVector neuron::fetch_EDF_parameters() const {
    NumericVector results = {A, tau, bias_term};
    results.names() = CharacterVector({"A", "tau", "bias_term"});
    return results;
  }

/*
 * ***********************************************************************************
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

double neuron::bounded_MSE_EDF_autocorr(
    const std::vector<double>& x, // 0 is A, 1 is tau
    std::vector<double>& grad,
    void* data                    // neuron object (this)
  ) {
    
    /*
     * EDF: autocor = A*exp(-lag/tau) + bias_term
     */
    
    // Grab neuron and estimated autocorrelation 
    neuron* nrn = static_cast<neuron*>(data);
    VectorXd est_autocorr = nrn->autocorr;
    int max_lag = est_autocorr.size();
    
    // Find mean squared error of autocorr prediction from these parameters x
    double bias_term = nrn->bias_term;
    double mse = 0.0;
    for (int i = 1; i < max_lag; i++) {
      // i = 1 because we throw out the first bin, which dominates
      double lag = (double)i;
      double pred = EDF_autocorr(lag, x[0], x[1], bias_term, 0);
      double err = est_autocorr(i) - pred;
      mse += err * err;
    }
    mse = mse/(double)(max_lag - 1.0);
    
    // Apply penalty to keep both parameters positive]
    double penalty_multiple = nrn->penalty_multiple;
    double p1 = penalty_multiple/(x[0] * x[0]);
    double p2 = penalty_multiple/(x[1] * x[1]);
    
    // Compute gradient if needed
    if (!grad.empty()) {
      
      // Compute gradient of mse
      std::vector<double> grd = {0.0, 0.0};
      for (int i = 1; i < max_lag; i++) {
        // i = 1 because we throw out the first bin, which dominates
        double lag = (double)i;
        double fx = est_autocorr(i) - EDF_autocorr(lag, x[0], x[1], bias_term, 0);
        double gA = -EDF_autocorr(lag, x[0], x[1], bias_term, 1);
        double gtau = -EDF_autocorr(lag, x[0], x[1], bias_term, 2);
        grd[0] += 2.0 * fx * gA;
        grd[1] += 2.0 * fx * gtau;
      }
      grd[0] = grd[0]/(double)(max_lag - 1.0);
      grd[1] = grd[1]/(double)(max_lag - 1.0);
      
      // Compute gradient of penalty
      grd[0] = grd[0] - (2.0 * penalty_multiple * x[0])/(x[0] * x[0] * x[0] * x[0]);
      grd[1] = grd[1] - (2.0 * penalty_multiple * x[1])/(x[1] * x[1] * x[1] * x[1]);
      
      grad.assign({grd[0], grd[1]});
    } 
    
    return mse + p1 + p2;
    
  }

void neuron::fit_autocorrelation() { 
   
    // Grab initial parameters
    std::vector<double> x = {A0, tau0};
    size_t n = x.size();
    
    // Set penalty multiple 
    int max_lag = autocorr.size();
    double sum_err0_sq = 0.0;
    for (int i = 1; i < max_lag; i++) {
      // i = 1 because we throw out the first bin, which dominates
      double lag = (double)i;
      double pred0 = EDF_autocorr(lag, A0, tau0, bias_term, 0);
      double err0 = autocorr(i) - pred0;
      sum_err0_sq += err0 * err0;
    }
    double mse0 = sum_err0_sq/(double)(max_lag - 1.0);
    
    // Set bias term
    bias_term = lambda * t_per_bin;
    bias_term = bias_term * bias_term;
    
    // When A is 1.0 and tau is 10.0, want p1 + p2 = mse0 * penalty_weight
    const double penalty_weight = 0.25;
    penalty_multiple = (mse0 * penalty_weight)/(1.0/(A0 * A0) + 1.0/(tau0 * tau0));
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LD_LBFGS, n); 
    opt.set_min_objective(neuron::bounded_MSE_EDF_autocorr, this);
    opt.set_ftol_rel(ctol);       // stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);   // Maximum number of evaluations to try
    
    // Fit model
    int success_code = 0;
    double min_fx;
    try {
      nlopt::result sc = opt.optimize(x, min_fx);
      success_code = static_cast<int>(sc);
    } catch (std::exception& e) { 
      if (false) {
        Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      } 
      success_code = 0;
    } 
    
    // Compute final modelled autocorrelation 
    autocorr_edf.resize(max_lag);
    for (int i = 0; i < max_lag; i++) {
      // Want i = 0 here for the construction of toeplitz sigma matrix
      double lag = (double)i;
      autocorr_edf(i) = EDF_autocorr(lag, x[0], x[1], bias_term, 0);
    }
    
    // Save optimization results, in proper time units
    A = x[0];
    tau = x[1] * t_per_bin;
    
  } 

double neuron::sigma_loss(
    const std::vector<double>& x, // the diagonal of the covariance matrix sigma
    std::vector<double>& grad,
    void* data                    // neuron object (this)
  ) {
    
    // Grab neuron 
    neuron* nrn = static_cast<neuron*>(data);
    
    // Grab fitted autocorrelation 
    NumericVector autocorr_fitted = to_NumVec(nrn->autocorr_edf);
    int max_lag = autocorr_fitted.size();
    
    // Grab threshold
    double dichot_threshold = nrn->gamma;
    
    // Extend sigma_diag and autocorr_fitted to have lag = 0 on front
    std::vector<double> sigma_diag(max_lag + 1);
    NumericVector autocorr_fitted0(max_lag + 1);
    for (int i = 0; i <= max_lag; i++) {
      if (i == 0) {
        sigma_diag[i] = 1.0;
        autocorr_fitted0(i) = 1.0;
      } else {
        sigma_diag[i] = x[i - 1];
        autocorr_fitted0(i) = autocorr_fitted(i - 1);
      }
    }
    
    // Create sigma from its diagonal
    NumericMatrix SIGMA = toeplitz(sigma_diag, sigma_diag);
    
    // Compute formula
    NumericVector output = dichot_gauss_sigma_formula(dichot_threshold, autocorr_fitted0, SIGMA);
    
    // Return sum of squares
    return Rcpp::sum(output * output);
    
  }

void neuron::dichot_gauss_parameters() {
    
    // Check that autocorrelation has been modeled
    if (autocorr_edf.size() == 0) {
      Rcpp::stop("autocorr_edf must be computed before dichot_gauss_parameters");
    }
    
    // Compute gamma (dichotomizing threshold) needed to simulate firing rate lambda 
    gamma = norm_cdf(
      1 - lambda,
      0,   // mean
      1,   // sd
      true // return inverse
    );
    
    // Use optimization to find the autocorrelation vector SIGMA needed to simulate the observed autocorrelation
    std::vector<double> x = to_dVec(autocorr_edf);
    size_t n = x.size();
    
    // Set up NLopt optimizer
    nlopt::opt opt(nlopt::LN_NELDERMEAD, n); // LD_LBFGS would need gradient function
    opt.set_min_objective(neuron::sigma_loss, this);
    opt.set_ftol_rel(ctol);       // stop when iteration changes objective fn value by less than this fraction 
    opt.set_maxeval(max_evals);   // Maximum number of evaluations to try
    
    // Find Sigma
    int success_code = 0;
    double min_fx;
    try {
      nlopt::result sc = opt.optimize(x, min_fx);
      success_code = static_cast<int>(sc);
    } catch (std::exception& e) {
      if (false) {
        Rcpp::Rcout << "Optimization failed: " << e.what() << std::endl;
      }
      success_code = 0;
    }
    
    // Save 
    sigma_gauss = x;
    
  }

neuron neuron::dichot_gauss_simulation(
    const int& trials
  ) {
    
    // Assuming that the binning was done by summation
    int max_lag = sigma_gauss.size();
    if (max_lag == 0) {Rcpp::stop("sigma_gauss must be computed before dichot_gauss_simulation");}
    
    // Add 1 to front of sigma_gauss 
    std::vector<double> sigma_gauss1(max_lag + 1);
    sigma_gauss1[0] = 1.0;
    for (int i = 0; i < max_lag; i++) {
      sigma_gauss1[i + 1] = sigma_gauss[i];
    }
    
    // Make random draws
    NumericMatrix sigma_gauss_matrix = toeplitz(sigma_gauss1, sigma_gauss1);
    NumericVector mu = rep(0.0, max_lag + 1);
    NumericMatrix simulated_trials_transpose = mvnorm_random(
      trials, 
      mu,
      sigma_gauss_matrix
    ); 
    
    // Take transpose, mvrnorm puts "points" (trials) as rows
    MatrixXd simulated_trials = to_eMat(simulated_trials_transpose).transpose();
    
    // Dichotomize
    for (int i = 0; i < simulated_trials.rows(); i++) {
      for (int j = 0; j < simulated_trials.cols(); j++) {
        simulated_trials(i, j) = (simulated_trials(i, j) < gamma) ? 0.0 : 1.0;
      }
    }
    
    // Make a copy of this neuron 
    neuron my_sim = neuron(*this);
    
    // Load with simulated trials
    // ... seems like these are needed for sim, but they also mess up the analysis of the observed data. 
    // my_sim.load_trial_data(simulated_trials);
    // my_sim.sim = true;
    // my_sim.unit_time = "bin";
    // my_sim.t_per_bin = 1.0;
    
    // Return simulated neuron
    return(my_sim);
    
  }

/*
 * RCPP_MODULE to expose class to R and function to initialize neuron
 */

RCPP_EXPOSED_CLASS(neuron)
RCPP_MODULE(neuron) {
  class_<neuron>("neuron")
  .constructor<int, std::string, std::string, std::string, bool, std::string, std::string, std::string, double, double>()
  .method("set_edf_initials", &neuron::set_edf_initials)
  .method("set_edf_termination", &neuron::set_edf_termination)
  .method("load_trial_data_R", &neuron::load_trial_data_R)
  .method("load_spike_raster_R", &neuron::load_spike_raster_R)
  .method("fetch_trial_data_R", &neuron::fetch_trial_data_R)
  .method("fetch_spike_raster_R", &neuron::fetch_spike_raster_R)
  .method("fetch_autocorr_R", &neuron::fetch_autocorr_R)
  .method("fetch_autocorr_edf_R", &neuron::fetch_autocorr_edf_R)
  .method("fetch_EDF_parameters", &neuron::fetch_EDF_parameters)
  .method("fetch_id_data", &neuron::fetch_id_data)
  .method("compute_autocorrelation", &neuron::compute_autocorrelation)
  .method("fit_autocorrelation", &neuron::fit_autocorrelation)
  .method("dichot_gauss_parameters", &neuron::dichot_gauss_parameters)
  .method("dichot_gauss_simulation", &neuron::dichot_gauss_simulation); // Returns MatrixXd??
}

