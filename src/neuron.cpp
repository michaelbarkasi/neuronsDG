
// neuron.cpp
#include "neuron.h"

/*
 * ***********************************************************************************
 * Helper functions
 */

const double Inf = 1e20;

// Build sequence of numbered string prefixes
CharacterVector enum_prefix(std::string prefix, int n) {
    CharacterVector result(n);
    for (int i = 0; i < n; ++i) {
      result[i] = prefix + std::to_string(i + 1);
    }
    return result;
  }

// Rolling mean
VectorXd roll_mean(
    const VectorXd& series,    // 1D vector of points to take rolling mean
    int filter_ws              // Size of window for taking rolling mean
  ) {
    int n = series.size();
    VectorXd series_out(n);
    for (int i = 0; i < n; i++) {
      int ws = std::min(i + 1, filter_ws);
      // Ensure iterator bounds stay within range
      auto start = series.begin() + (i - ws + 1);
      auto end = series.begin() + (i + 1);
      series_out[i] = std::accumulate(start, end, 0.0) / static_cast<double>(ws);
    } 
    return series_out;
  }
// ... overload
NumericVector roll_mean(
    const NumericVector& series,    
    int filter_ws  
  ) {
    int n = series.size();
    NumericVector series_out(n);
    for (int i = 0; i < n; i++) {
      int ws = std::min(i + 1, filter_ws);
      // Ensure iterator bounds stay within range
      auto start = series.begin() + (i - ws + 1);
      auto end = series.begin() + (i + 1);
      series_out[i] = std::accumulate(start, end, 0.0) / static_cast<double>(ws);
    } 
    return series_out;
  }
// ... overload
std::vector<double> roll_mean(
    const std::vector<double>& series,
    int filter_ws   
  ) {
    int n = series.size();
  std::vector<double> series_out(n);
    for (int i = 0; i < n; i++) {
      int ws = std::min(i + 1, filter_ws);
      // Ensure iterator bounds stay within range
      auto start = series.begin() + (i - ws + 1);
      auto end = series.begin() + (i + 1);
      series_out[i] = std::accumulate(start, end, 0.0) / static_cast<double>(ws);
    } 
    return series_out;
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
// ... overload 
NumericVector to_NumVec(
    const std::vector<double>& vec
  ) {
    return wrap(vec);
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

// Empirical Pearson correlation between two vectors
double empirical_corr(
    const VectorXd& x,
    const VectorXd& y
  ) {
   
    // Check dimensions
    int n = x.size();
    if (y.size() != n) {Rcpp::stop("Vectors must be same length for empirical correlation calculation");}
    
    // Calculate components
    double mean_x = x.sum()/(double)n;
    double mean_y = y.sum()/(double)n;
    double mean_xy = x.dot(y)/(double)n;
    double mean_x2 = x.dot(x)/(double)n;
    double mean_y2 = y.dot(y)/(double)n;
    double sd_x = sqrt(mean_x2 - mean_x * mean_x);
    double sd_y = sqrt(mean_y2 - mean_y * mean_y);
    
    // Return correlation estimate
    double corr = (mean_xy - mean_x * mean_y)/(sd_x * sd_y);
    return corr;
    
  }

// Empirical Pearson correlation between two variables sampled many times 
double empirical_corr_multisample(
    const MatrixXd& X, // Rows as intratrial samples, columns as trials
    const MatrixXd& Y  // Rows as intratrial samples, columns as trials
  ) {
    
    // Check dimensions
    int n_trials = X.cols();
    if (Y.cols() != n_trials) {Rcpp::stop("Matrices must have same number of columns (trials) for empirical correlation calculation");}
    
    // Calculate correlation within each trial
    double corr = 0.0;
    int n_trials_nonan = n_trials;
    for (int n = 0; n < n_trials; n++) {
      double corr_n = empirical_corr(X.col(n), Y.col(n));
      if (std::isnan(corr_n)) {
        n_trials_nonan -= 1;
      } else {
        corr += corr_n;
      }
    }
    
    // Return mean correlation across trials
    return corr/(double)n_trials_nonan;
    
  }

// Estimate Pearson correlation across lags
VectorXd empirical_corr_lagged(
    const MatrixXd& TS1, // Time series 1, rows as time points, columns as trials
    const MatrixXd& TS2  // Time series 2, rows as time points, columns as trials
  ) {
    
    // Check for matching dimensions 
    const int max_lag = TS1.rows();
    if (TS2.rows() != max_lag) {Rcpp::stop("Time series must have same number of time points (rows)");}
    const int N_trial = TS1.cols();
    if (TS2.cols() != N_trial) {Rcpp::stop("Time series must have same number of trials (columns)");}
    
    // Initialize vector to hold autocorrelation values
    VectorXd lagged_corr(max_lag);
    lagged_corr.setZero();
    
    // Compute correlation for all possible lags
    for (int lag = 0; lag < max_lag; lag++) {
      MatrixXd TS2_shifted = TS2(seq(0 + lag, max_lag), Eigen::all);
      MatrixXd TS1_cut = TS1(seq(0, max_lag - lag), Eigen::all);
      lagged_corr(lag) = empirical_corr_multisample(TS1_cut, TS2_shifted);
    }
    
    // Return correlation 
    return(lagged_corr); 
    
  }

// Estimate raw correlation across lags, raw version (no mean subtraction, no normalization by std)
VectorXd empirical_corr_lagged_raw(
    const MatrixXd& TS1, // Time series 1, rows as time points, columns as trials
    const MatrixXd& TS2  // Time series 2, rows as time points, columns as trials
  ) {
    
    // Check for matching dimensions 
    const int max_lag = TS1.rows();
    if (TS2.rows() != max_lag) {Rcpp::stop("Time series must have same number of time points (rows)");}
    const int N_trial = TS1.cols();
    if (TS2.cols() != N_trial) {Rcpp::stop("Time series must have same number of trials (columns)");}
    
    // Initialize vector to hold autocorrelation values
    VectorXd lagged_corr(max_lag);
    lagged_corr.setZero();
    
    // Compute correlation for all possible lags
    for (int lag = 0; lag < max_lag; lag++) {
      
      // Shift and cut
      MatrixXd TS2_shifted = TS2(seq(0 + lag, max_lag), Eigen::all);
      MatrixXd TS1_cut = TS1(seq(0, max_lag - lag), Eigen::all);
      
      // Sum over trials
      for (int n = 0; n < N_trial; n++) {
        lagged_corr(lag) += TS2_shifted.col(n).dot(TS1_cut.col(n));
      }
      // The normalization term is the number of overlapping samples at this lag. It's more transparent 
      //   to define it in terms of max_bin instead of max_lag, but in this case max_bin == max_lag.
      double normalization_term = 1.0/((double)max_lag - (double)lag);
      lagged_corr(lag) *= normalization_term;
      // Find mean proportion of co-active bins per trial
      lagged_corr(lag) *= (1.0/N_trial);
      
    }
    
    // Return correlation 
    return(lagged_corr); 
    
  }

// EDF model 
double EDF_autocorr(
    const double& lag,
    const double& A, 
    const double& tau, 
    const double& bias_term,
    const int& return_grad
  ) {
    if (return_grad == 0) {        // No gradient, return function output
      return A * exp(-lag/tau) + bias_term;
    } else if (return_grad == 1) { // gradient wrt A
      return exp(-lag/tau);
    } else if (return_grad == 2) { // gradient wrt tau
      return A * exp(-lag/tau) * (lag/(tau*tau));
    } else {
      return 0.0;
    }
  }

// Multivariate normal CDF, upper tail
double mvnorm_cdf_uppertail(
    const NumericVector& threshold, 
    const NumericMatrix& sigma  // covariance matrix of dimension n less than 1000
  ) {
    
    if (sigma.nrow() != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (sigma.nrow() >= 1000) {Rcpp::stop("Covariance matrix must be less than 1000x1000");}
    if (sigma.nrow() != threshold.size()) {Rcpp::stop("Matrix diagonal and threshold vector must be same length");}
    
    Function pmvnorm("pmvnorm", Environment::namespace_env("mvtnorm"));
    // ... uses Genz algorithm
    
    double prob = as<double>(
      pmvnorm( // by default, lower = -Inf, upper = Inf, and mean = 0.
        Named("lower") = threshold, 
        Named("sigma") = sigma, 
        Named("keepAttr") = false
      )
    );
    
    return prob;
    
  }

// Normal CDF, with inverse
double norm_cdf(
    const double& x,
    const double& mu,
    const double& sd,
    const bool& inverse
  ) {
    double xc = x;
    using boost::math::normal; 
    normal standard_normal(mu, sd);
    if (inverse) {
      if (xc < 1e-10) {xc = 1e-10;}
      if (xc > 1 - 1e-10) {xc = 1.0 - 1e-10;}
      return boost::math::quantile(standard_normal, xc);
    } else {
      return boost::math::cdf(standard_normal, xc);
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
NumericVector dg_sigma_formula(
    const double& threshold,      // threshold for dichotomization
    const NumericVector& cov,     // desired covarance after dichotomization
    const NumericMatrix& sigma    // covariance matrix of multivariate Gaussian
  ) {
    // We know threshold and cov. By finding the sigma which sends this function 
    //   to zero, we can find the covariance needed for dichotomized Gaussian simulation
   
    // Check dimension
    int dim = sigma.nrow();
    if (dim != sigma.ncol()) {Rcpp::stop("Covariance matrix must be square");}
    if (dim != cov.size()) {Rcpp::stop("Covariance vector must have the same length as sigma diagonal");}
    
    // Find probability of a point being above the threshold along all dimensions
    NumericVector threshold_vec = Rcpp::rep(threshold, dim); 
    double Phi2_upper = mvnorm_cdf_uppertail(
      threshold_vec, 
      sigma
    );
   
    // Find probability of a point being below the threshold along one dimension
    double Phi = norm_cdf(
      threshold, 
      0.0,     // mean
      1.0,     // sd
      false    // return inverse? No, return cdf
    );
    
    // desired sigma will be the one which sends all elements to zero
    //  Formula: cov = Phi2_upper - (1 - Phi) * (1 - Phi) is derived as follows: 
    //   By definition of cov, cov = E[X1*X2] - E[X1]*E[X2].
    //   In this case, E[X1] = E[X2] = P(X > threshold) = 1 - Phi.
    //   X1*X2 != 0 only if both X1 and X2 > threshold, which occurs with probability Phi2_upper.
    NumericVector residuals(dim);
    for (int i = 0; i < dim; i++) {
      residuals[i] = cov[i] - Phi2_upper + (1 - Phi) * (1 - Phi);
    }
    
    return residuals;
    
  }

// Wrapper for use with find-root-bisection algorithm 
double dg_sigma_formula_scalar(
    const double& threshold,      // threshold for dichotomization
    const double& cov,            // desired covarance after dichotomization
    const double& sigma           // Gaussian covariance
  ) {
    
    // Construct covariance matrix sigma (2x2)
    NumericMatrix sigmaMat(2, 2);
    sigmaMat(_,0) = NumericVector::create(1.0, sigma);
    sigmaMat(_,1) = NumericVector::create(sigma, 1.0);
    
    // Include self-covariance on front, which is the variance, sd^2
    //   The Gaussian is normal, so mean is zero and sd is 1, so variance is 1
    NumericVector cov1 = {1.0, cov};
    
    // Evaluate formula and return second value
    NumericVector residual = dg_sigma_formula(threshold, cov1, sigmaMat);
    // Return only the second element, corresponding to cov
    return residual[1]; 
    
  }

// Function to find sigma by root bisection 
double dg_find_sigma_RootBisection(
    const double& threshold,      // threshold for dichotomization
    const double& cov             // desired covarance after dichotomization
  ) {
    
    // Set search parameters 
    const int max_iter = 50; 
    const double tol = 1e-4;
    
    // Initiate sigmas
    double sigma_lower = -0.999;
    double sigma_upper = 0.999;
    
    // Evaluate formula
    double fx_lower = dg_sigma_formula_scalar(threshold, cov, sigma_lower);
    double fx_upper = dg_sigma_formula_scalar(threshold, cov, sigma_upper);
    
    // Run checks 
    if (abs(fx_lower) < tol) {return sigma_lower;}
    else if (abs(fx_upper) < tol) {return sigma_upper;}
    else if (fx_lower * fx_upper > tol) {return 0.0;} // Both initial covariance values lie on same side of zero crossing
    
    // Run bisection
    double fx = Inf;
    double sigma_mid;
    int iter = 0;
    while (abs(fx) > tol && iter < max_iter) {
      
      // Find midpoint
      sigma_mid = (sigma_lower + sigma_upper)/2.0;
      fx = dg_sigma_formula_scalar(threshold, cov, sigma_mid);
      
      // Update bounds
      if (fx > 0.0) {
        sigma_lower = sigma_mid;
      } else {
        sigma_upper = sigma_mid;
      }
      
      // Update iteration
      iter++;
      
    }
    
    return sigma_mid; 
    
  }

// Function to make a matrix positive definite
NumericMatrix makePositiveDefinite(
    const NumericMatrix& NumX
  ) {
    
    MatrixXd X = to_eMat(NumX);
    SelfAdjointEigenSolver<MatrixXd> solver(X);
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();
    
    // Ensure all eigenvalues are positive
    for (int i = 0; i < eigenvalues.size(); ++i) {
      if (eigenvalues(i) < 1e-10) {  // Adjust small or negative values
        eigenvalues(i) = 1e-10;
      }
    }
    
    // Reconstruct the matrix
    return to_NumMat(eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose());
    
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
    const std::string genotype,
    const std::string sex,
    const std::string hemi,
    const std::string region,
    const std::string age,
    const bool sim, 
    const std::string unit_time, 
    const std::string unit_sample_rate, 
    const std::string unit_data, 
    const double t_per_bin, 
    const double sample_rate
  ) : id_num(id_num), 
    recording_name(recording_name), 
    type(type), 
    genotype(genotype),
    sex(sex),
    hemi(hemi), 
    region(region),
    age(age),
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
 * Member function implementations, basic data handling, loading
 */

// Loading trial data (Eigen matrix)
void neuron::load_trial_data(const MatrixXd& td) {
    
    // Set trial data
    trial_data = td;
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
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
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
    // Make spike raster
    if (unit_data == "spike") {
      infer_raster();
    }
    
  } 

// Loading spike raster (Eigen matrix)
void neuron::load_spike_raster(const MatrixXd& sr) {
    
    // Check input form
    if (unit_data != "spike") {Rcpp::stop("Spike raster can only be loaded if unit_data is 'spike'");}
    if (sr.cols() != 2) {Rcpp::stop("Spike raster must have two columns");}
   
    // Save spike raster
    spike_raster = sr;
    
    // Infer trial data
    infer_trial();
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
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
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
  } 

/*
 * ***********************************************************************************
 * Member function implementations, basic data handling, fetching
 */

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

// Method to return fields with neuron ID information
List neuron::fetch_id_data() const {
    
    return List::create(
      _["id_num"] = id_num,
      _["recording_name"] = recording_name,
      _["type"] = type,
      _["genotype"] = genotype,
      _["sex"] = sex,
      _["hemi"] = hemi,
      _["region"] = region,
      _["age"] = age,
      _["sim"] = sim,
      _["unit_time"] = unit_time,
      _["unit_sample_rate"] = unit_sample_rate,
      _["unit_data"] = unit_data,
      _["t_per_bin"] = t_per_bin,
      _["sample_rate"] = sample_rate
    );
    
  }

// Method to return firing rate
NumericVector neuron::fetch_lambda() const {
    NumericVector lambdas = {lambda, lambda_bin};
    lambdas.names() = CharacterVector({"lambda", "lambda_bin"});
    return lambdas;
  }

/*
 * ***********************************************************************************
 * Member function implementations, basic data handling, misc
 */

// Infer trial data from spike raster
void neuron::infer_trial() {
    
    // Find needed dimensions
    int num_trials = static_cast<int>(spike_raster.col(1).maxCoeff());
    int max_time = static_cast<int>(spike_raster.col(0).maxCoeff() - spike_raster.col(0).minCoeff());
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

/*
 * ***********************************************************************************
 * Member function implementations, fetch analysis results
 */

// Fetching autocorrelation 
VectorXd neuron::fetch_autocorr() const {return autocorr;}

// Fetching autocorrelation 
NumericVector neuron::fetch_autocorr_R() const {return wrap(autocorr);}
NumericVector neuron::fetch_autocorr_edf_R() const {return wrap(autocorr_edf);}
NumericVector neuron::fetch_sigma_gauss_R() const {return wrap(sigma_gauss);}

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
    const std::string& bin_count_action,
    const bool& use_raw
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
      
      // Compute standard deviation
      spike_sd_bin = sqrt(data.array().square().sum()/(double)(data.rows() * data.cols()) - lambda_bin * lambda_bin);
      
    } else {
      data = trial_data;
      spike_sd_bin = spike_sd;
    }
    
    // Find autocorrelation
    if (use_raw) {
      autocorr = empirical_corr_lagged_raw(data, data);
    } else {
      autocorr = empirical_corr_lagged(data, data);
    }
    
    autocorr = roll_mean(autocorr, t_per_bin); // smooth with rolling mean
    
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
    
    // Apply penalty to keep both parameters positive
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
    nlopt::srand(static_cast<unsigned long>(bias_term*1e6));
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

void neuron::dg_parameters(
    const bool& use_raw,
    const bool& verbose
  ) {
    
    if (verbose) {Rcpp::Rcout << "Finding dichotomized Gaussian parameters for neuron " << id_num << ", " << recording_name << " ..." << std::endl;}
    
    // Check that autocorrelation has been modeled
    int max_lag = autocorr_edf.size();
    if (max_lag == 0) {
      Rcpp::stop("autocorr_edf must be computed before dg_parameters");
    }
    
    // Compute gamma (dichotomizing threshold) needed to simulate firing rate lambda 
    gamma = norm_cdf(
      1 - lambda_bin,
      0.0,   // mean
      1.0,   // sd
      true   // return inverse
    );
    
    // Resize sigma_gauss
    sigma_gauss.resize(max_lag);
    
    // Find the covariance sigma_gauss for each lag
    if (use_raw) {
      for (int i = 0; i < max_lag; i++) {
        // When raw, autocorr_edf = E[X1*X2], so need to subtract off E[X1]*E[X2] = lambda^2 to get covariance
        // Want lambda_bin (not lambda) because these parameters are going to simulations, which have time units of bin
        sigma_gauss[i] = dg_find_sigma_RootBisection(gamma, autocorr_edf[i] - lambda_bin*lambda_bin);
      }
    } else {
      for (int i = 0; i < max_lag; i++) {
        // When using Pearson correlation, autocorr_edf = cov/(sd1*sd2), so need to multiply by sd1*sd2 to get covariance
        // Want spike_sd_bin (not lambda) because these parameters are going to simulations, which have time units of bin
        sigma_gauss[i] = dg_find_sigma_RootBisection(gamma, autocorr_edf[i] * spike_sd_bin*spike_sd_bin);
      }
    }
    
   
  }

// Make dichotomized Gaussian simulation of neuron
neuron neuron::dg_simulation(
    const int& trials,
    const bool& verbose
  ) {
   
    if (verbose) {Rcpp::Rcout << "Running dichotomized Gaussian simulation of neuron " << id_num << ", " << recording_name << " ..." << std::endl;}
    
    // Check that the covariance matrix sigma_guass has been computed
    int max_lag = sigma_gauss.size();
    if (max_lag == 0) {Rcpp::stop("sigma_gauss must be computed before dg_simulation");}
    
    // Make random draws
    NumericMatrix sigma_gauss_matrix = makePositiveDefinite(toeplitz(sigma_gauss, sigma_gauss));
    NumericVector mu = rep(0.0, max_lag); // Drawing from MVN with mean = 0 and sd = 1
    NumericMatrix simulated_trials_transpose = mvnorm_random(
      trials, 
      mu,
      sigma_gauss_matrix // Must be covariance matrix (positive definite)
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
    my_sim.load_trial_data(simulated_trials);
    my_sim.sim = true;
    my_sim.unit_time = "bin";
    my_sim.t_per_bin = 1.0;
    
    // return simulation 
    return my_sim;
    
  }

// Estimate autocorrelation parameters for a neuron with dichotomized Gaussian simulations
NumericMatrix neuron::estimate_autocorr_params(
    const int& trials_per_sim, 
    const int& num_sims,
    const std::string& bin_count_action,
    const double& A0,
    const double& tau0,
    const double& ctol,
    const int& max_evals,
    const bool& use_raw,
    const bool& verbose
  ) {
    
    if (verbose) {
      Rcpp::Rcout << "Estimating autocorrelation parameters for neuron " << id_num << ", " << recording_name << " ..." << std::endl;
    }
    
    // Pre-process base neuron, if needed
    if (autocorr_edf.size() == 0) {
      // ... set exponential decay function (EDF) parameters 
      set_edf_initials(A0, tau0);
      set_edf_termination(ctol, max_evals);
      // ... compute autocorrelation
      compute_autocorrelation(bin_count_action, use_raw);
      // ... fit autocorrelation 
      fit_autocorrelation();
    }
    
    // Find parameters for dichotomized Gaussian simulation
    dg_parameters(
      use_raw,
      false // verbose?
    );
    
    // Initialize matrix to hold results 
    NumericMatrix sim_results(num_sims, 9); 
    
    for (int s = 0; s < num_sims; s++) {
      
      // Simulate neuron with dichotomized Gaussian
      neuron my_sim = dg_simulation(
        trials_per_sim, 
        false // verbose?
      );
      
      // Set exponential decay function (EDF) parameters 
      my_sim.set_edf_initials(A0, tau0);
      my_sim.set_edf_termination(ctol, max_evals);
      
      // Compute autocorrelation
      my_sim.compute_autocorrelation(bin_count_action, use_raw);
      
      // Fit autocorrelation 
      my_sim.fit_autocorrelation();
      
      // Fetch results 
      VectorXd autocorr_edf_tail_eigen = my_sim.autocorr_edf.tail(my_sim.autocorr_edf.size() - 1);
      NumericVector autocorr_edf_tail = to_NumVec(autocorr_edf_tail_eigen);
      NumericVector sim_results_row = {
        my_sim.lambda, 
        my_sim.lambda_bin, 
        my_sim.A, 
        my_sim.tau * t_per_bin, // convert back to original time units
        my_sim.bias_term,
        my_sim.autocorr_edf[0],
        max(autocorr_edf_tail),
        mean(autocorr_edf_tail),
        min(autocorr_edf_tail)
      };
      
      sim_results.row(s) = sim_results_row;
      
    }
    
    return sim_results;
    
  }

/*
 * RCPP_MODULE to expose class to R and function to initialize neuron
 */

RCPP_EXPOSED_CLASS(neuron)
RCPP_MODULE(neuron) {
  class_<neuron>("neuron")
  .constructor<int, std::string, std::string, std::string, std::string, std::string, std::string, std::string, bool, std::string, std::string, std::string, double, double>()
  .method("set_edf_initials", &neuron::set_edf_initials)
  .method("set_edf_termination", &neuron::set_edf_termination)
  .method("load_trial_data_R", &neuron::load_trial_data_R)
  .method("load_spike_raster_R", &neuron::load_spike_raster_R)
  .method("fetch_trial_data_R", &neuron::fetch_trial_data_R)
  .method("fetch_spike_raster_R", &neuron::fetch_spike_raster_R)
  .method("fetch_autocorr_R", &neuron::fetch_autocorr_R)
  .method("fetch_autocorr_edf_R", &neuron::fetch_autocorr_edf_R)
  .method("fetch_sigma_gauss_R", &neuron::fetch_sigma_gauss_R)
  .method("fetch_EDF_parameters", &neuron::fetch_EDF_parameters)
  .method("fetch_id_data", &neuron::fetch_id_data)
  .method("fetch_lambda", &neuron::fetch_lambda)
  .method("compute_autocorrelation", &neuron::compute_autocorrelation)
  .method("fit_autocorrelation", &neuron::fit_autocorrelation)
  .method("dg_parameters", &neuron::dg_parameters)
  .method("dg_simulation", &neuron::dg_simulation)
  .method("estimate_autocorr_params", &neuron::estimate_autocorr_params);
}

