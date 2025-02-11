
// neuron.h
#ifndef NEURON_H
#define NEURON_H

// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <nlopt.hpp>
using namespace Rcpp;
using namespace Eigen;

// Helper functions
CharacterVector enum_prefix(std::string prefix, int n);

// Neuron class

class neuron {

  private:
    
    // ID parameters
    const int id_num = 0;                               // Fixed ID number for each neuron
    const std::string recording_name = "not_provided";  // Recording (if any) on which this neuron is based
    const std::string type = "generic";                 // Type of neuron, e.g. "generic", "blackbox" "LIF", "McCullochPitts", "excitatory", "inhibitory", etc.
    const std::string hemi = "not_provided";            // Hemisphere of neuron, e.g. "left", "right"
    const bool sim = false;                             // Whether this neuron is simulated or based on recorded data
    
    // Unit specifications
    const std::string unit_time = "ms";                 // Unit of time, e.g., "ms", "bin", "sample"
    const std::string unit_sample_rate = "Hz";          // Unit of recording sample rate, e.g., "Hz", "kHz"
    const std::string unit_data = "mV";                 // Unit of data, e.g., "mV", "spike count"
    
    // Unit conversions 
    const double t_per_bin = 1.0;                       // Time (in above units) per bin, e.g., 1 ms per bin
    const double sample_rate = 1e4;                     // Sample rate (in above units), e.g., 10000 Hz
    
    // Data fields
    MatrixXd trial_data;                                // NxM matrix of doubles, rows as recording times (in "unit_time"), columns as trials, data values in "unit_data"
    MatrixXd spike_raster;                              // Nx2 matrix, each row one spike, columns as time (in "unit_time") and trial number
    double lambda;                                      // Mean value of neuron, in "unit_data" per "unit_time"
    
    // Analysis fields
    VectorXd autocorr;                                  // Estimated (observed) autocorrelation of trial_data
    VectorXd autocorr_edf;                              // Estimated autocorrelation of trial_data using EDF model
    double A0 = 1.0;                                    // Initial amplitude of EDF model of autocorrelation 
    double tau0 = 10.0;                                 // Initial time constant of EDF model of autocorrelation
    double A;                                           // Fitted amplitude of EDF model of autocorrelation
    double tau;                                         // Fitted time constant of EDF model of autocorrelation
    double bias_term;                                   // Bias term for EDF model of autocorrelation
    double penalty_multiple;                            // For scaling boundary penalty terms when fitting EDF model of autocorrelation

  public:
    
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
    
    // Member functions for fetching analysis results
    VectorXd fetch_autocorr() const;
    NumericVector fetch_autocorr_R() const;
    NumericVector fetch_autocorr_edf_R() const;
    
    // Member functions for data analysis
    void compute_autocorrelation(const std::string& bin_count_action); // action must be 'sum', 'boolean', or 'mean'
    static double bounded_MSE_EDF_autocorr(
        const std::vector<double>& x, // 0 is A, 1 is tau
        std::vector<double>& grad,
        void* data                    // neuron object (this)
    );
    void fit_autocorrelation();

};

#endif
