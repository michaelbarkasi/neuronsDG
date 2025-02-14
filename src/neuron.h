
// neuron.h
#ifndef NEURON_H
#define NEURON_H

// Rcpp
// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <nlopt.hpp>
#include <boost/math/distributions/normal.hpp>
using namespace Rcpp;
using namespace Eigen;

// Helper functions

// Build sequence of numbered string prefixes
CharacterVector enum_prefix(std::string prefix, int n);

// Conert between vector types
std::vector<double> to_dVec(const VectorXd& vec);
VectorXd to_eVec(const std::vector<double>& vec);

// Exponential decay function (with gradients) for modelling
double EDF_autocorr(
    const double& lag, 
    const double& A, 
    const double& tau, 
    const double& bias_term, 
    const int& return_grad // 0 = function output, 1 = gradient wrt A, 2 = gradient wrt tau
  );

// Probability distributions
double mvnorm_cdf(
    const VectorXd& x, 
    const VectorXd& mu, 
    const MatrixXd& sigma
  );

double norm_cdf(
    const double& x, 
    const double& mu, 
    const double& sd,
    const bool& inverse
  );

MatrixXd mvnorm_random(
    int n, 
    VectorXd mu, 
    MatrixXd sigma
  );

// Function to create Toeplitz matrix
NumericMatrix toeplitz(
    const std::vector<double>& first_col, 
    const std::vector<double>& first_row
  );

// Formula for estimating sigma for dichotomized Gaussian simulation
NumericVector dichot_gauss_sigma_formula(
    const double& threshold,      // threshold for dichotomization
    const NumericVector& cov,     // desired covarance after dichotomization
    const NumericMatrix& sigma    // covariance matrix
  );

// Neuron class

class neuron {
  
  // private:

  public:
    
    // ID parameters
    int id_num = 0;                               // Fixed ID number for each neuron
    std::string recording_name = "not_provided";  // Recording (if any) on which this neuron is based
    std::string type = "generic";                 // Type of neuron, e.g. "generic", "blackbox" "LIF", "McCullochPitts", "excitatory", "inhibitory", etc.
    std::string hemi = "not_provided";            // Hemisphere of neuron, e.g. "left", "right"
    bool sim = false;                             // Whether this neuron is simulated or based on recorded data
    
    // Unit specifications
    std::string unit_time = "ms";                 // Unit of time, e.g., "ms", "bin", "sample"
    std::string unit_sample_rate = "Hz";          // Unit of recording sample rate, e.g., "Hz", "kHz"
    std::string unit_data = "mV";                 // Unit of data, e.g., "mV", "spike count"
    
    // Unit conversions 
    double t_per_bin = 1.0;                       // Time (in above units) per bin, e.g., 1 ms per bin
    double sample_rate = 1e4;                     // Sample rate (in above units), e.g., 10000 Hz
    
    // Data fields
    MatrixXd trial_data;                          // NxM matrix of doubles, rows as recording times (in "unit_time"), columns as trials, data values in "unit_data"
    MatrixXd spike_raster;                        // Nx2 matrix, each row one spike, columns as time (in "unit_time") and trial number
    double lambda;                                // Mean value of neuron, in "unit_data" per "unit_time"
    double lambda_bin;                            // Mean value of neuron, in "unit_data" per bin
    
    // Analysis fields
    VectorXd autocorr;                            // Estimated (observed) autocorrelation of trial_data
    VectorXd autocorr_edf;                        // Estimated autocorrelation of trial_data using EDF model
    std::vector<double> sigma_gauss;              // Covariance matrix for Gaussian simulation
    int max_evals = 500;                          // Max number of evals when fitting EDF to autocorrelation
    double ctol = 1e-7;                           // Convergence tolerance when fitting EDF to autocorrelation
    double A0 = 0.1;                              // Initial amplitude of EDF model of autocorrelation 
    double tau0 = 10.0;                           // Initial time constant of EDF model of autocorrelation
    double A;                                     // Fitted amplitude of EDF model of autocorrelation
    double tau;                                   // Fitted time constant of EDF model of autocorrelation
    double bias_term;                             // Bias term for EDF model of autocorrelation
    double penalty_multiple;                      // For scaling boundary penalty terms when fitting EDF model of autocorrelation
    double gamma;                                 // Threshold for dichotomized Gaussian simulation
    
    // Constructor and Destructor
    neuron(
      const int id_num = 0, 
      const std::string recording_name = "not_provided", 
      const std::string type = "generic", 
      const std::string hemi = "not_provided",
      const bool sim = false, 
      const std::string unit_time = "ms", 
      const std::string unit_sample_rate = "Hz", 
      const std::string unit_data = "mV", 
      const double t_per_bin = 1.0, 
      const double sample_rate = 1e4
    );
    virtual ~neuron() {};
    
    // Member functions for adjusting settings
    void set_edf_initials(double a0, double t0);
    void set_edf_termination(double ct, int me);
    
    // Member functions for loading data
    void load_trial_data(const MatrixXd& td);
    void load_trial_data_R(const NumericMatrix& td);
    void load_spike_raster(const MatrixXd& sr);
    void load_spike_raster_R(const NumericMatrix& sr);
    void infer_raster();
    void infer_trial();
    
    // Member functions for fetching data
    MatrixXd fetch_trial_data() const;
    NumericMatrix fetch_trial_data_R() const;
    MatrixXd fetch_spike_raster() const;
    NumericMatrix fetch_spike_raster_R() const;
    List fetch_id_data() const;
    NumericVector fetch_lambda() const;
    
    // Member functions for fetching analysis results
    VectorXd fetch_autocorr() const;
    NumericVector fetch_autocorr_R() const;
    NumericVector fetch_autocorr_edf_R() const;
    NumericVector fetch_sigma_gauss_R() const; // continuous (normal) equivalent of autocorr_edf
    NumericVector fetch_EDF_parameters() const;
    
    // Member functions for data analysis
    void compute_autocorrelation(
      const std::string& bin_count_action // action must be 'sum', 'boolean', or 'mean'
    ); 
    static double bounded_MSE_EDF_autocorr(
      // Objective function for fitting EDF model to autocorrelation
      const std::vector<double>& x, // 0 is A, 1 is tau
      std::vector<double>& grad,
      void* data                    // neuron object (this)
    );
    void fit_autocorrelation();
    static double sigma_loss(
      const std::vector<double>& x,
      std::vector<double>&grad,
      void* data
    );
    void dichot_gauss_parameters();
    neuron dichot_gauss_simulation(const int& trials);

};

#endif
