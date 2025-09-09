
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

// Convert between vector types
std::vector<double> to_dVec(const VectorXd& vec);
std::vector<double> to_dVec(const NumericVector& vec);
VectorXd to_eVec(const std::vector<double>& vec);
NumericVector to_NumVec(const VectorXd& vec);
NumericVector to_NumVec(const std::vector<double>& vec);

// Exponential decay function (with gradients) for modelling
double EDF_autocorr(
    const double& lag, 
    const double& A, 
    const double& tau, 
    const double& bias_term, 
    const int& return_grad // 0 = function output, 1 = gradient wrt A, 2 = gradient wrt tau
  );

// multivariate normal CDF
double mvnorm_cdf(
    const NumericVector& upper, 
    const NumericMatrix& sigma
  );

// Normal CDF inverse
double norm_cdf(
    const double& x, 
    const double& mu, 
    const double& sd,
    const bool& inverse
  );

// Multivariate normal random number generator
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
NumericVector dg_sigma_formula(
    const double& threshold,      // threshold for dichotomization
    const NumericVector& cov,     // desired covarance after dichotomization
    const NumericMatrix& sigma    // covariance matrix
  );

// Wrapper for use with find-root-bisection algorithm 
double dg_sigma_formula_scalar(
    const double& threshold,      // threshold for dichotomization
    const double& cov,            // desired covarance after dichotomization
    const double& sigma           // covariance matrix
  );

// Function to find sigma by root bisection 
double dg_find_sigma_RootBisection(
    const double& threshold,      // threshold for dichotomization
    const double& cov             // desired covarance after dichotomization
  );

// Function to make a matrix positive definite
NumericMatrix makePositiveDefinite(
    const NumericMatrix& NumX
  );

// Neuron class

class neuron {
  
  // private:
  
  // public:

  public:
    
    // ID parameters
    int id_num = 0;                               // Fixed ID number for neuron
    std::string recording_name = "not_provided";  // Recording (if any) on which this neuron is based
    std::string type = "generic";                 // Modeled electrophysiology of neuron, e.g. "generic", "blackbox" "LIF", "McCullochPitts", "excitatory", "inhibitory", etc.
    std::string genotype = "WT";                  // Genotype of animal, e.g. "WT", "KO", "MECP2", "transgenic", etc.
    std::string sex = "not_provided";             // Sex of animal
    std::string hemi = "not_provided";            // Hemisphere of neuron, e.g. "left", "right"
    std::string region = "not_provided";          // Brain region of neuron, e.g. "V1", "M1", "CA1", "PFC", etc.
    std::string age = "not_provided";             // Age of animal, e.g. "P0", "P7", "P14", "adult", etc.
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
    std::vector<double> sigma_gauss;              // Covariance values for Gaussian simulation, in time units of bin
    int max_evals = 500;                          // Max number of evals when fitting EDF to autocorrelation
    double ctol = 1e-7;                           // Convergence tolerance when fitting EDF to autocorrelation
    double A0 = 0.1;                              // Initial amplitude of EDF model of autocorrelation 
    double tau0 = 10.0;                           // Initial time constant of EDF model of autocorrelation
    double A;                                     // Fitted amplitude of EDF model of autocorrelation
    double tau;                                   // Fitted time constant of EDF model of autocorrelation
    double bias_term;                             // Bias term for EDF model of autocorrelation
    double penalty_multiple;                      // For scaling boundary penalty terms when fitting EDF model of autocorrelation
    double gamma;                                 // Threshold for dichotomized Gaussian simulation (in time units of bin)
    
    // Constructor and Destructor
    neuron(
      const int id_num = 0, 
      const std::string recording_name = "not_provided", 
      const std::string type = "generic", 
      const std::string genotype = "WT",
      const std::string sex = "not_provided",
      const std::string hemi = "not_provided",
      const std::string region = "not_provided",
      const std::string age = "not_provided",
      const bool sim = false, 
      const std::string unit_time = "ms", 
      const std::string unit_sample_rate = "Hz", 
      const std::string unit_data = "mV", 
      const double t_per_bin = 1.0, 
      const double sample_rate = 1e4
    );
    virtual ~neuron() {};
    
    // Copy method 
    neuron(const neuron& other) = default;
    
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
    void dg_parameters(const bool& verbose);
    neuron dg_simulation(const int& trials, const bool& verbose);
    NumericMatrix estimate_autocorr_params(
        const int& trials_per_sim, 
        const int& num_sims,
        const std::string& bin_count_action,
        const double& A0,
        const double& tau0,
        const double& ctol,
        const int& max_evals,
        const bool& verbose
    );

};

#endif
