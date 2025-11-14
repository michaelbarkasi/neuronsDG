
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
    const VectorXd& series,         // 1D vector of points to take rolling mean
    int filter_ws                   // Size of window for taking rolling mean
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

// Return logical vector giving elements of left which match right
LogicalVector eq_left_broadcast(
    const CharacterVector& left,
    const String& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }
// ... overload
LogicalVector eq_left_broadcast(
    const std::vector<int>& left,
    const int& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }
// ... overload
LogicalVector eq_left_broadcast(
    const VectorXi& left,
    const int& right
  ) {
    int n = left.size();
    LogicalVector out(n);
    for (int i = 0; i < n; i++) {
      out[i] = left[i] == right;
    }
    return out;
  }

// Convert boolean masks to integer indexes
IntegerVector Rwhich(
    const LogicalVector& x
  ) {
    std::vector<int> indices;  // Use std::vector for efficient dynamic resizing
    for (int i = 0; i < x.size(); ++i) {
      if (x[i]) {
        indices.push_back(i);
      }
    }
    return wrap(indices);  // Convert std::vector to IntegerVector
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

// Convert to Eigen matrix with integers
MatrixXi to_eiMat(
    const IntegerMatrix& X
  ) {
    int Xnrow = X.nrow();
    int Xncol = X.ncol();
    MatrixXi M = MatrixXi(Xnrow, Xncol);
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
// ... overload
NumericMatrix to_NumMat(
    const MatrixXi& M
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

// Convert to IntegerMatrix
IntegerMatrix to_IntMat(
    const MatrixXi& M
  ) {
    int M_nrow = M.rows();
    int M_ncol = M.cols();
    IntegerMatrix X(M_nrow, M_ncol);
    for (int i = 0; i < M_nrow; i++) {
      for (int j = 0; j < M_ncol; j++) {
        X(i, j) = M(i, j);
      }
    }
    return X;
  }

/*
 * ***********************************************************************************
 * Functions for computing correlations
 */

// Empirical Pearson correlation between two vectors
double empirical_corr(
    const VectorXd& x,
    const VectorXd& y,
    const bool& use_raw
  ) {
   
    // Check dimensions
    int n = x.size();
    if (y.size() != n) {Rcpp::stop("Vectors must be same length for empirical correlation calculation");}
    
    // Calculate components for raw correlation
    double mean_xy = x.dot(y)/(double)n;
    if (use_raw) {
      return mean_xy;
    } else {
      
      // Calculate components for Pearson correlation
      double mean_x = x.sum()/(double)n;
      double mean_y = y.sum()/(double)n;
      double mean_x2 = x.dot(x)/(double)n;
      double mean_y2 = y.dot(y)/(double)n;
      double sd_x = sqrt(mean_x2 - mean_x * mean_x);
      double sd_y = sqrt(mean_y2 - mean_y * mean_y);
      
      // Return correlation estimate
      double corr = (mean_xy - mean_x * mean_y)/(sd_x * sd_y);
      return corr;
      
    }
    
  }

// Empirical Pearson correlation between two variables sampled many times 
double empirical_corr_multisample(
    const MatrixXd& X,            // Rows as intratrial samples, columns as trials
    const MatrixXd& Y,            // Rows as intratrial samples, columns as trials
    const bool& use_raw
  ) {
    
    // Check dimensions
    int n_trials = X.cols();
    if (Y.cols() != n_trials) {Rcpp::stop("Matrices must have same number of columns (trials) for empirical correlation calculation");}
    // Calculate correlation within each trial
    double corr = 0.0;
    int n_trials_nonan = n_trials;
    for (int n = 0; n < n_trials; n++) {
      double corr_n = empirical_corr(X.col(n), Y.col(n), use_raw);
      if (std::isnan(corr_n) || std::isinf(corr_n)) {
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
    const MatrixXd& TS_ref,       // Reference time series, rows as time points, columns as trials
    const MatrixXd& TS_cmp,       // Comparison time series (to be lagged), rows as time points, columns as trials
    const int& max_lag,
    const bool& use_raw
  ) {
    
    // Check for matching dimensions 
    const int n_ref = TS_ref.rows();
    const int n_cmp = TS_cmp.rows();
    const int N_trial = TS_ref.cols();
    if (max_lag > n_ref) {Rcpp::stop("max_lag cannot be greater than number of time points (rows) in the reference time series");}
    if (n_cmp < n_ref) {Rcpp::stop("Comparison time series must have at least as many rows as the reference series");}
    if (TS_cmp.cols() != N_trial) {Rcpp::stop("Time series must have same number of trials (columns)");}
    // Check if padding by zeros is needed
    MatrixXd TS_cmp_padded;
    if (n_cmp < n_ref + max_lag) {
      TS_cmp_padded = MatrixXd::Zero(n_ref + max_lag, N_trial);
      TS_cmp_padded(seq(0, n_cmp - 1), Eigen::all) = TS_cmp;
    } else {
      TS_cmp_padded = TS_cmp; 
    }
    
    // Initialize vector to hold autocorrelation values
    VectorXd lagged_corr(max_lag);
    lagged_corr.setZero();
    
    // Compute correlation for all possible lags
    for (int lag = 0; lag < max_lag; lag++) {
      MatrixXd TS_cmp_shifted = TS_cmp_padded(seq(0 + lag, n_ref + lag - 1), Eigen::all);
      lagged_corr(lag) = empirical_corr_multisample(TS_ref, TS_cmp_shifted, use_raw);
    }
    
    // Return correlation 
    return(lagged_corr); 
    
  }

// Estimate raw correlation across lags, raw version (no mean subtraction, no normalization by std)
VectorXd empirical_corr_lagged_raw(
    const MatrixXd& TS1,          // Time series 1, rows as time points, columns as trials
    const MatrixXd& TS2           // Time series 2, rows as time points, columns as trials
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

/*
 * ***********************************************************************************
 * Probability distribution functions
 */

// Multivariate normal CDF, upper tail
double mvnorm_cdf_uppertail(
    const NumericVector& threshold, 
    const NumericMatrix& sigma    // covariance matrix of dimension n less than 1000
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
    int n,                        // Number of points to generate
    NumericVector mu,             // Mean vector, length determines dimension
    NumericMatrix sigma           // Covariance matrix, square, same dimension as mu
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

/*
 * ***********************************************************************************
 * Dichotomized Gaussian helper functions
 */

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

/*
 * ***********************************************************************************
 * Matrix and vector operations
 */

// Create Toeplitz matrix
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
    return to_NumMat(MatrixXd(eigenvectors * eigenvalues.asDiagonal() * eigenvectors.transpose()));
    
  }

/*
 * ***********************************************************************************
 * Neuron and network classes
 */

// Constructor, neuron
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

// Constructor, motif
motif::motif(
    const std::string motif_name
  ) : motif_name(motif_name)
  { 
      // No initialization operations
  }

// Constructor, network
network::network(
    const std::string network_name, 
    const std::string recording_name, 
    const std::string type, 
    const std::string genotype,
    const std::string sex,
    const std::string hemi,
    const std::string region,
    const std::string age,
    const std::string unit_time, 
    const std::string unit_sample_rate, 
    const std::string unit_potential, 
    const std::string unit_current,
    const std::string unit_conductance,
    const std::string unit_distance,
    const double t_per_bin, 
    const double sample_rate
  ) : network_name(network_name), 
    recording_name(recording_name), 
    type(type), 
    genotype(genotype),
    sex(sex),
    hemi(hemi), 
    region(region),
    age(age),
    unit_time(unit_time), 
    unit_sample_rate(unit_sample_rate), 
    unit_potential(unit_potential), 
    unit_current(unit_current),
    unit_conductance(unit_conductance),
    unit_distance(unit_distance),
    t_per_bin(t_per_bin), 
    sample_rate(sample_rate)
  { 
      // No initialization operations
  }

/*
 * ***********************************************************************************
 * Neuron and network member function implementations, adjust settings
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

void network::set_network_structure(
    CharacterVector nrn_types,
    std::vector<int> nrn_type_valences,
    CharacterVector lyr_names,
    int n_lyr,
    int n_cls,
    double lyr_height,
    double cls_width,
    double lyr_separation_factor,
    double cls_separation_factor,
    MatrixXi nrn_per_node,
    std::vector<MatrixXd> recur_factors,
    double pruning_thresh_factor
  ) {
    
    // Check layer names (needed for motifs)
    if (lyr_names.size() != n_lyr) {
      Rcpp::Rcout << "lyr_names size: " << lyr_names.size() << ", n_layers: " << n_lyr << std::endl;
      Rcpp::stop("Length of lyr_names must equal n_layers");
    }
   
    // Set parameters
    neuron_types = nrn_types;
    neuron_type_valences = nrn_type_valences;
    layer_names = lyr_names;
    n_layers = n_lyr;
    n_columns = n_cls;
    layer_height = lyr_height;
    column_width = cls_width;
    layer_separation_factor = lyr_separation_factor;
    column_separation_factor = cls_separation_factor;
    neurons_per_node = nrn_per_node;
    recurrence_factors = recur_factors;
    pruning_threshold_factor = pruning_thresh_factor;
    
    // Set network components
    n_neuron_types = neuron_types.size();
    int n_nodes = n_layers * n_columns;
    n_neurons = 0; // Compute total number of neurons as we go
    node_range_ends.assign(n_nodes, 0);
    node_coordinates_spatial.resize(n_nodes, 2);
    for (int l = 0; l < n_layers; l++) {
      for (int c = 0; c < n_columns; c++) {
        int node_idx = l * n_columns + c;
        // Set global spatial coordinates for this node
        node_coordinates_spatial(node_idx, 0) = c * column_width * column_separation_factor;
        node_coordinates_spatial(node_idx, 1) = l * layer_height * layer_separation_factor;
        for (int t = 0; t < n_neuron_types; t++) {
          // Randomly select neuron numbers for each node
          int n = (int)R::rpois(neurons_per_node(l,t));
          // Keep track of the number of cells assigned so far
          n_neurons += n; 
          // Keep track of the types of these cells
          for (int i = 0; i < n; i++) {
            neuron_type_name.push_back(neuron_types[t]);
            neuron_type_num.push_back(t);
          }
        }
        // Save end-point index for this node
        node_range_ends[node_idx] = n_neurons - 1;
      }
    }
    
    // Resize network coordinate components 
    coordinates_spatial = MatrixXd::Zero(n_neurons, 2); 
    coordinates_node = MatrixXi::Zero(n_neurons, 2); // layer, column
    
  };

void network::set_network_structure_R(
    CharacterVector neuron_types,
    std::vector<int> neuron_type_valences,
    CharacterVector layer_names,
    int n_layers,
    int n_columns,
    double layer_height,
    double column_width,
    double layer_separation_factor,
    double column_separation_factor,
    IntegerMatrix neurons_per_node,
    List recurrence_factors,
    double pruning_threshold_factor
  ) {
    
    // Convert recurrence factors from R List to std::vector<MatrixXd>
    std::vector<MatrixXd> rec_factors_vec;
    for (int i = 0; i < recurrence_factors.size(); i++) {
      NumericMatrix rec_mat_r = recurrence_factors[i];
      MatrixXd rec_mat = to_eMat(rec_mat_r);
      rec_factors_vec.push_back(rec_mat);
    }
    
    // Call main function
    set_network_structure(
      neuron_types,
      neuron_type_valences,
      layer_names,
      n_layers,
      n_columns,
      layer_height,
      column_width,
      layer_separation_factor,
      column_separation_factor,
      to_eiMat(neurons_per_node),
      rec_factors_vec,
      pruning_threshold_factor
    );
    
  };

/*
 * ***********************************************************************************
 * Network member function implementations, build network
 */

void motif::load_projection(
    const Projection& proj,
    const int& max_up,
    const int& max_down,
    const double& c_strength
  ) {
    projections.push_back(proj);
    max_col_shift_up.push_back(max_up);
    max_col_shift_down.push_back(max_down);
    connection_strength.push_back(c_strength);
    n_projections++;
  }

// Function to set transconductances and spatial coordinates for all local nodes 
void network::make_local_nodes() {
    
    if (edge_types.size() != 0) {
      Rcpp::Rcout << "Edge types have already been set; cannot run make_local_nodes twice; returning." << std::endl;
      return;
    }
    
    // Initialize vectors to track local edge coordinates
    std::vector<int> local_edges_pre; 
    std::vector<int> local_edges_post;
    
    // Initialize local transconductance matrix
    MatrixXd local_transconductances = MatrixXd::Zero(n_neurons, n_neurons);
   
    // Layer index of the local node
    for (int l = 0; l < n_layers; l++) {
      
      // Get recurrence factor matrix for this layer
      MatrixXd recurrence_factor_matrix = recurrence_factors[l];
      
      // Column index of the local node
      for (int c = 0; c < n_columns; c++) {
        
        // Get node ID number
        int node_idx = l * n_columns + c;
        // Get spatial position of this node
        double node_x = node_coordinates_spatial(node_idx, 0);
        double node_y = node_coordinates_spatial(node_idx, 1);
        // Get the range of neuron ID numbers for this node
        int node_range_start = (node_idx == 0) ? 0 : node_range_ends[node_idx - 1] + 1;
        int node_range_end = node_range_ends[node_idx];
        
        // For all combinations of pre- and post-synaptic neurons in this node
        for (int idx_pre = node_range_start; idx_pre <= node_range_end; idx_pre++) {
          
          // Set spatial coordinates
          coordinates_spatial(idx_pre, 0) = node_x + R::rnorm(0.0, column_width);
          coordinates_spatial(idx_pre, 1) = node_y + R::rnorm(0.0, layer_height);
          
          // Set node coordinates
          coordinates_node(idx_pre, 0) = c;
          coordinates_node(idx_pre, 1) = l;
          
          // Set transconductance into post-synaptic cells
          for (int idx_post = node_range_start; idx_post <= node_range_end; idx_post++) {
            
            // Get neuron types for pre- and post-synaptic neurons
            int t_pre = neuron_type_num[idx_pre];
            int t_post = neuron_type_num[idx_post];
            
            // Get neuron valences for pre-synaptic neurons
            double val_pre = (double)neuron_type_valences[t_pre];
            
            // Get recurrence factor for this connection type
            double rec_factor = recurrence_factor_matrix(t_post, t_pre);
            double pruning_threshold = pruning_threshold_factor * rec_factor;
            
            // Set transductance 
            double trans = std::abs(R::rnorm(0.0, rec_factor));
            if (trans > pruning_threshold) {
              local_transconductances(idx_post, idx_pre) = val_pre * trans;
              // Save edge coordinate
              local_edges_pre.push_back(idx_pre);
              local_edges_post.push_back(idx_post);
            }
            
          }
          
        }
        
      }
    }
    
    // Save to transconductance matrix
    transconductances.push_back(local_transconductances);
    
    // Collect local edge coordinates in matrix
    int n_local_edges = local_edges_pre.size();
    MatrixXi local_edges(n_local_edges, 2); 
    local_edges.col(0) = Eigen::Map<VectorXi>(local_edges_pre.data(), n_local_edges);
    local_edges.col(1) = Eigen::Map<VectorXi>(local_edges_post.data(), n_local_edges);
    
    // Save to edge types
    edge_types.push_back(local_edges);
    
  }

// Function to apply circuit motif
void network::apply_circuit_motif(
    const motif& cmot
  ) {
   
    if (edge_types.size() != 1) {
      Rcpp::Rcout << "Must set local edges before applying any circuit motifs; returning." << std::endl;
      return;
    }
    
    // Initialize vectors to track motif edge coordinates
    std::vector<int> motif_edges_pre; 
    std::vector<int> motif_edges_post;
    
    // Initialize motif transconductance matrix
    MatrixXd motif_transconductances = MatrixXd::Zero(n_neurons, n_neurons);
    
    // For each projection in the motif
    for (int p = 0; p < cmot.n_projections; p++) {
      
      // Grab projection
      Projection proj = cmot.projections[p];
      
      // Grab projection strength and set pruning threshold
      double proj_strength = cmot.connection_strength[p];
      double pruning_threshold = pruning_threshold_factor * proj_strength;
      
      // Grab pre and post cell types 
      int type_pre = Rwhich(eq_left_broadcast(neuron_types, proj.pre_type))[0];
      int type_post = Rwhich(eq_left_broadcast(neuron_types, proj.post_type))[0];
      // ... and make masks 
      LogicalVector pre_type_mask = eq_left_broadcast(neuron_type_num, type_pre);
      LogicalVector post_type_mask = eq_left_broadcast(neuron_type_num, type_post);
      
      // Grab pre-synaptic valence
      int val_pre = neuron_type_valences[type_pre];
      
      // Grab pre and post layers
      int layer_pre = Rwhich(eq_left_broadcast(layer_names, proj.pre_layer))[0];
      int layer_post = Rwhich(eq_left_broadcast(layer_names, proj.post_layer))[0];
      // ... and make masks 
      LogicalVector pre_layer_mask = eq_left_broadcast(coordinates_node(Eigen::all,1), layer_pre);
      LogicalVector post_layer_mask = eq_left_broadcast(coordinates_node(Eigen::all,1), layer_post);
      
      // Grab pre and post densities
      double density_pre = proj.pre_density;
      double density_post = proj.post_density;
      
      // Grab max shifts
      int max_up = cmot.max_col_shift_up[p];
      int max_down = cmot.max_col_shift_down[p];
      
      // Build column range
      VectorXi col_range(max_down - max_down + 1);
      for (int i = 0; i < col_range.size(); i++) col_range[i] = max_down + i;
      
      // Pre-make all column masks 
      LogicalMatrix column_masks(n_neurons, n_columns);
      for (int c = 0; c < n_columns; c++) {
        column_masks(_, c) = eq_left_broadcast(coordinates_node(Eigen::all,0), c);
      }
      
      // Apply projection to each column 
      for (int c = 0; c < n_columns; c++) {
        
        // Get pre-synaptic column mask
        LogicalVector pre_column_mask = column_masks(_, c);
        
        // Shift range to this column 
        VectorXi col_range_shifted = col_range.array() + c;
        
        // For each target column
        for (int tc : col_range_shifted) {
          
          // Check if target column is valid
          if (tc < 0 || tc >= n_columns) {continue;}
          
          // Get post-synaptic column mask
          LogicalVector post_column_mask = column_masks(_, tc);
          
          // Sample pre-synaptic cells 
          LogicalVector pre_mask = pre_type_mask & pre_layer_mask & pre_column_mask;
          IntegerVector pre_indices = Rwhich(pre_mask);
          int n_pre = R::rpois(pre_indices.size() * density_pre);
          n_pre = std::min(n_pre, 1);
          IntegerVector pre_sampled = Rcpp::sample(pre_indices, n_pre, false);
          
          // Prepare for repeated sampling of post-synaptic cells
          LogicalVector post_mask = post_type_mask & post_layer_mask & post_column_mask;
          IntegerVector post_indices = Rwhich(post_mask);
          
          for (int pre_c : pre_sampled) {
            
            // Sample post-synaptic cells
            int n_post = R::rpois(post_indices.size() * density_post);
            n_post = std::min(n_post, 1);
            IntegerVector post_sampled = Rcpp::sample(post_indices, n_post, false);
            
            // Set transductances 
            for (int post_c : post_sampled) {
              double trans = std::abs(R::rnorm(0.0, proj_strength));
              if (trans > pruning_threshold) {
                motif_transconductances(post_c, pre_c) = val_pre * trans;
                // Save edge coordinate
                motif_edges_pre.push_back(pre_c);
                motif_edges_post.push_back(post_c);
              }
              
            }
            
          }
          
        }
        
      }
      
    }
    
    // Save to transconductance matrix vector 
    transconductances.push_back(motif_transconductances);
    
    // Collect local edge coordinates in matrix
    int n_motif_edges = motif_edges_pre.size();
    MatrixXi motif_edges(n_motif_edges, 2); 
    motif_edges.col(0) = Eigen::Map<VectorXi>(motif_edges_pre.data(), n_motif_edges);
    motif_edges.col(1) = Eigen::Map<VectorXi>(motif_edges_post.data(), n_motif_edges);
    
    // Save to edge types
    edge_types.push_back(motif_edges);
    
    // Add motif name
    edge_type_names.push_back(cmot.motif_name);
    
  }

/*
 * ***********************************************************************************
 * Neuron member function implementations, basic data handling, loading
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
void neuron::load_spike_raster(
    const MatrixXd& sr,
    const int& min_duration,
    const int& max_displacement
  ) {
    
    // Check input form
    if (unit_data != "spike") {Rcpp::stop("Spike raster can only be loaded if unit_data is 'spike'");}
    if (sr.cols() != 2) {Rcpp::stop("Spike raster must have two columns");}
   
    // Save spike raster
    spike_raster = sr;
    
    // Infer trial data
    infer_trial(min_duration, max_displacement);
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
  } 

// Loading spike raster (Rcpp matrix)
void neuron::load_spike_raster_R(
    const NumericMatrix& sr,
    const int& min_duration,
    const int& max_displacement
  ) {
    
    // Check input form
    if (unit_data != "spike") {Rcpp::stop("Spike raster can only be loaded if unit_data is 'spike'");}
    if (sr.ncol() != 2) {Rcpp::stop("Spike raster must have two columns");}
    CharacterVector sr_colnames = colnames(sr); 
    if (sr_colnames[0] != "time" || sr_colnames[1] != "trial") {Rcpp::stop("Spike raster must have columns named 'time' and 'trial'");}
    
    // Save spike raster
    spike_raster = as<MatrixXd>(sr);
    
    // Infer trial data
    infer_trial(min_duration, max_displacement);
    
    // Compute mean neuron value (e.g., firing rate)
    lambda = trial_data.sum()/(double)(trial_data.rows() * trial_data.cols());
    lambda_bin = lambda * t_per_bin;
    
    // Compute standard deviation
    spike_sd = sqrt(trial_data.array().square().sum()/(double)(trial_data.rows() * trial_data.cols()) - lambda * lambda);
    
  } 

/*
 * ***********************************************************************************
 * Neuron and network member function implementations, basic data handling, fetching
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

// Method to fetch network components 
List network::fetch_network_components() const {
   
    // Convert transconductances into list of NumericMatrix
    List transconductance_matrices(transconductances.size());
    for (int tci = 0; tci < transconductances.size(); tci++) {
      MatrixXd tc = transconductances[tci];
      NumericMatrix tc_r = to_NumMat(tc);
      transconductance_matrices[tci] = tc_r;
    } 
    
    // Convert edge_types into list of NumericMatrix
    List edge_type_matrices(edge_types.size());
    CharacterVector emn = CharacterVector::create("pre_neuron_idx", "post_neuron_idx");
    for (int eti = 0; eti < edge_types.size(); eti++) {
      MatrixXi et = edge_types[eti];
      NumericMatrix et_r = to_NumMat(et);
      colnames(et_r) = emn;
      edge_type_matrices[eti] = et_r;
    }
    
    return List::create(
      _["n_neurons"] = n_neurons,
      _["n_neuron_types"] = n_neuron_types,
      _["transconductances"] = transconductance_matrices,
      _["node_coordinates_spatial"] = to_NumMat(node_coordinates_spatial),
      _["coordinates_spatial"] = to_NumMat(coordinates_spatial),
      _["coordinates_node"] = to_NumMat(coordinates_node),
      _["neuron_type_name"] = neuron_type_name,
      _["neuron_type_num"] = neuron_type_num,
      _["node_range_ends"] = node_range_ends,
      _["edge_types"] = edge_type_matrices, 
      _["edge_type_names"] = edge_type_names
    );
   
  }

/*
 * ***********************************************************************************
 * Neuron member function implementations, basic data handling, misc
 */

// Infer trial data from spike raster
void neuron::infer_trial(
    const int& min_duration,
    const int& max_displacement
  ) {
    
    // Find number of trials
    int num_trials = static_cast<int>(spike_raster.col(1).maxCoeff());
    
    // Find duration
    int duration = static_cast<int>(spike_raster.col(0).maxCoeff() - spike_raster.col(0).minCoeff());
    duration = ((duration + 9) / 10) * 10; // round up to nearest 10
    if (duration < min_duration) {duration = min_duration;}
    if (duration <= 0) {Rcpp::stop("Spike raster has zero or negative duration");}
    
    // Find displacement
    int displacement = static_cast<int>(std::round(spike_raster.col(0).minCoeff()));
    if (displacement > max_displacement) {displacement = max_displacement;}
    
    // Resize trial_data and set to zero
    trial_data.resize(duration, num_trials);
    trial_data.setZero();
    
    // Fill trial data with recorded spikes 
    for (int r = 0; r < spike_raster.rows(); r++) {
      int time_idx = std::max(0, static_cast<int>(std::round(spike_raster(r, 0))) - displacement - 1);
      int trial_idx = static_cast<int>(std::round(spike_raster(r, 1))) - 1;
      trial_data(time_idx, trial_idx) += 1.0;
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
 * Neuron member function implementations, fetch analysis results
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
 * Neuron member function implementations, analysis
 */

// Compute cross-correlation of this neuron with another neuron
VectorXd neuron::compute_crosscorrelation(
    const neuron& nrn_compare,
    const std::string& bin_count_action,
    const int& max_lag,
    const bool& use_raw,
    const bool& verbose
  ) {
    
    if (verbose) {
      Rcpp::Rcout << "Computing cross-correlation between neuron " << id_num << " (ref) and neuron " << nrn_compare.id_num << " (comparison)" << std::endl;
    }
    
    // Compare time units 
    if (unit_time != nrn_compare.unit_time) {
      Rcpp::stop("Both neurons must have the same unit_time for cross-correlation calculation");
    } else if (verbose) {
      Rcpp::Rcout << "Both neurons have unit_time = " << unit_time << std::endl;
    }
    
    // Compare bin size 
    if (t_per_bin != nrn_compare.t_per_bin) {
      Rcpp::stop("Both neurons must have the same t_per_bin for cross-correlation calculation");
    } else if (verbose) {
      Rcpp::Rcout << "Both neurons have t_per_bin = " << t_per_bin << " " << unit_time << std::endl;
    }
    
    // Compare trial lengths
    const int trial_length_ref = trial_data.rows();
    const int trial_length_comp = nrn_compare.trial_data.rows();
    if (verbose) {
      Rcpp::Rcout << "Reference neuron trial length: " << trial_length_ref << " " << unit_time << std::endl;
      Rcpp::Rcout << "Comparison neuron trial length: " << trial_length_comp << " " << unit_time << std::endl;
    }
    if (trial_length_comp < trial_length_ref) {Rcpp::stop("Comparison neuron must have trial length at least as long as reference neuron for cross-correlation calculation");}
    if (max_lag > trial_length_ref) {Rcpp::stop("max_lag must be less than or equal to the trial length of the reference neuron for cross-correlation calculation");}
    
    // Check time unit of the trial matrix, home cell
    bool convert_to_bins = false;
    int max_lag_bin = max_lag;
    int trial_length_ref_bin = trial_length_ref;
    int trial_length_comp_bin = trial_length_comp;
    if (unit_time != "bin") {
      convert_to_bins = true;
      max_lag_bin = (int) max_lag/t_per_bin;
      trial_length_ref_bin = (int) trial_length_ref/t_per_bin;
      trial_length_comp_bin = (int) trial_length_comp/t_per_bin;
      if (verbose) {
        Rcpp::Rcout << "Converting to bins for cross-correlation calculation" << std::endl;
        Rcpp::Rcout << "max_lag in bins: " << max_lag_bin << std::endl;
        Rcpp::Rcout << "Reference neuron trial length in bins: " << trial_length_ref_bin << std::endl;
        Rcpp::Rcout << "Comparison neuron trial length in bins: " << trial_length_comp_bin << std::endl;
      }
    }
    
    // Check binning action
    if (bin_count_action != "sum" && bin_count_action != "boolean" && bin_count_action != "mean") {
      Rcpp::stop("bin_count_action must be 'sum', 'boolean', or 'mean'");
    } else if (verbose) {
      Rcpp::Rcout << "Using bin_count_action = " << bin_count_action << std::endl;
    }
    
    // Initialize matrices for data manipulation 
    const int N_trial = std::min(trial_data.cols(), nrn_compare.trial_data.cols());
    MatrixXd binned_data_ref(trial_length_ref_bin, N_trial);
    MatrixXd binned_data_comp(trial_length_comp_bin, N_trial);
    
    // Collapse bins if needed
    if (convert_to_bins) {
      
      for (int n = 0; n < N_trial; n++) {
        for (int b = 0; b < trial_length_ref_bin; b++) {
          // Sum spikes in bin
          double bin_count_ref = trial_data(seq(b*(int)t_per_bin, (b*(int)t_per_bin) + t_per_bin - 1), n).sum();
          // Apply binning action
          if (bin_count_action == "mean") {
            binned_data_ref(b, n) = bin_count_ref/(double)t_per_bin;
          } else if (bin_count_action == "sum") {
            binned_data_ref(b, n) = bin_count_ref;
          } else if (bin_count_action == "boolean") {
            binned_data_ref(b, n) = (bin_count_ref > 0.0) ? 1.0 : 0.0;
          }
        }
        for (int b = 0; b < trial_length_comp_bin; b++) {
          // Sum spikes in bin
          double bin_count_comp = nrn_compare.trial_data(seq(b*(int)t_per_bin, (b*(int)t_per_bin) + t_per_bin - 1), n).sum();
          // Apply binning action
          if (bin_count_action == "mean") {
            binned_data_comp(b, n) = bin_count_comp/(double)t_per_bin;
          } else if (bin_count_action == "sum") {
            binned_data_comp(b, n) = bin_count_comp;
          } else if (bin_count_action == "boolean") { 
            binned_data_comp(b, n) = (bin_count_comp > 0.0) ? 1.0 : 0.0;
          }
        } 
      }
      
    } else {
      // No binning
      binned_data_ref = trial_data;
      binned_data_comp = nrn_compare.trial_data;
    }
    
    // How many times should the lag be stepped across the data?
    int n_lag_steps = trial_length_ref_bin - max_lag_bin;
    int n_lag_steps_good = n_lag_steps;
    
    // compute cross-correlation as average of lagged correlations across all possible lag steps
    if (verbose) {
      Rcpp::Rcout << "Computing cross-correlation vector of length " << max_lag_bin << " bins" << std::endl;
    }
    VectorXd crosscorr(max_lag_bin);
    crosscorr.setZero();
    for (int lag_step = 0; lag_step <= n_lag_steps; lag_step++) {
      
      // Subset data to current lag step
      int ref_end = lag_step + max_lag_bin - 1;
      MatrixXd binned_data_ref_cut = binned_data_ref(seq(lag_step, ref_end), Eigen::all);
      // Comparison data should be longer, as it is shifted by max_lag_bin
      int comp_end = std::min(ref_end + max_lag_bin, trial_length_comp_bin);
      MatrixXd binned_data_comp_cut = binned_data_comp(seq(lag_step, comp_end), Eigen::all);
      // Find lagged correlation for this step
      VectorXd step_crosscorr = empirical_corr_lagged(binned_data_ref_cut, binned_data_comp_cut, max_lag_bin, use_raw);
      
      // If lagged correlation is well defined, add to running sum
      if ((step_crosscorr.array().isNaN()).any() || (step_crosscorr.array().isInf()).any()) {
        n_lag_steps_good -= 1;
      } else {
        crosscorr += step_crosscorr;
      }
    }
    // Reduce running sum to average
    crosscorr /= (double)(n_lag_steps + 1);
    
    if (verbose) {
      Rcpp::Rcout << "Possible lag steps: " << n_lag_steps << ", good lag steps: " << n_lag_steps_good << std::endl;
      Rcpp::Rcout << "Cross-correlation calculation complete" << std::endl;
    }
    
    crosscorr = roll_mean(crosscorr, t_per_bin); // smooth with rolling mean
    
    // Return 
    return crosscorr;
    
  }

// Wrapper
NumericVector neuron::compute_crosscorrelation_R(
    const neuron& nrn_compare,
    const std::string& bin_count_action,
    const int& max_lag,
    const bool& use_raw,
    const bool& verbose
  ) {
    VectorXd cc = compute_crosscorrelation(nrn_compare, bin_count_action, max_lag, use_raw, verbose);
    return to_NumVec(cc);
  }

// Compute autocorrelation of trial data
void neuron::compute_autocorrelation(
    const std::string& bin_count_action,
    int max_lag,
    const bool& use_raw
  ) {
    
    /*
     * Based on Neophytou et al 2022, "Differences in temporal processing ... "
     *   https://doi.org/10.1371/journal.pbio.3001803
     *   equation: x_n(t) = trial_matrix[t,n]
     */
    
    // Check time unit of the trial matrix
    bool convert_to_bins = false;
    int T_n = trial_data.rows();
    if (max_lag > T_n) {
      std::cout << "T_n: " << T_n << std::endl;
      std::cout << "max_lag: " << max_lag << std::endl;
      Rcpp::stop("max_lag must be less than or equal to the trial length");
    }
    if (unit_time != "bin") {
      convert_to_bins = true;
      max_lag = (int) max_lag/t_per_bin;
      T_n = (int) T_n/t_per_bin;
    }
    
    // Check binning action
    if (bin_count_action != "sum" && bin_count_action != "boolean" && bin_count_action != "mean") {
      Rcpp::stop("bin_count_action must be 'sum', 'boolean', or 'mean'");
    }
    
    // Initialize matrix for data manipulation 
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
      
      // Compute standard deviation, used below, in dg_parameters
      spike_sd_bin = sqrt(data.array().square().sum()/(double)(data.rows() * data.cols()) - lambda_bin * lambda_bin);
      
    } else {
      data = trial_data;
      spike_sd_bin = spike_sd;
    }
    
    // Find autocorrelation
    // if (use_raw) {
    //   autocorr = empirical_corr_lagged_raw(data, data);
    // } else {
    //   autocorr = empirical_corr_lagged(data, data, use_raw);
    // }
    autocorr = empirical_corr_lagged(data, data, max_lag, use_raw);
    
    autocorr = roll_mean(autocorr, t_per_bin); // smooth with rolling mean
    
  }

// Objective function for fitting EDF model to autocorrelation
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

// Fit EDF model to autocorrelation
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

// Find parameters for dichotomized Gaussian simulation
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
    int max_lag,
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
      compute_autocorrelation(bin_count_action, max_lag, use_raw);
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
      my_sim.compute_autocorrelation(bin_count_action, max_lag, use_raw);
      
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
 * RCPP_MODULE to expose class to R
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
  .method("compute_crosscorrelation_R", &neuron::compute_crosscorrelation_R)
  .method("fit_autocorrelation", &neuron::fit_autocorrelation)
  .method("dg_parameters", &neuron::dg_parameters)
  .method("dg_simulation", &neuron::dg_simulation)
  .method("estimate_autocorr_params", &neuron::estimate_autocorr_params);
}

RCPP_EXPOSED_CLASS(motif)
RCPP_MODULE(motif) {
  class_<motif>("motif")
  .constructor<std::string>()
  .method("load_projection", &motif::load_projection);
}

RCPP_EXPOSED_CLASS(network)
RCPP_MODULE(network) {
  class_<network>("network")
  .constructor<std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, std::string, double, double>()
  .method("set_network_structure", &network::set_network_structure)
  .method("set_network_structure_R", &network::set_network_structure_R)
  .method("make_local_nodes", &network::make_local_nodes)
  .method("apply_circuit_motif", &network::apply_circuit_motif)
  .method("fetch_network_components", &network::fetch_network_components);
}

RCPP_EXPOSED_CLASS(Projection)
RCPP_MODULE(Projection) {
  class_<Projection>("Projection")
  .constructor()
  .field("pre_type",      &Projection::pre_type)
  .field("pre_layer",     &Projection::pre_layer)
  .field("pre_density",   &Projection::pre_density)
  .field("post_type",     &Projection::post_type)
  .field("post_layer",    &Projection::post_layer)
  .field("post_density",  &Projection::post_density);
}