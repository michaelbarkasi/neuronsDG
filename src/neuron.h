
// neuron.h
#ifndef NEURON_H
#define NEURON_H

// Rcpp
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
using namespace Eigen;

class neuron {

  private:
    
    // ID parameters
    const int id_num = 0;                               // Fixed ID number for each neuron
    const std::string recording_name = "not_provided";  // Recording (if any) on which this neuron is based
    const std::string type = "generic";                 // Type of neuron, e.g. "generic", "blackbox" "LIF", "McCullochPitts", "excitatory", "inhibitory", etc.
    const bool sim = false;                             // Whether this neuron is simulated or based on recorded data
    
    // Unit specifications
    const std::string unit_time = "ms";                 // Unit of time, e.g., "ms", "bin", "sample"
    const std::string unit_sample_rate = "Hz";          // Unit of recording sample rate, e.g., "Hz", "kHz"
    const std::string unit_data = "mV";                 // Unit of data, e.g., "mV", "spike count"
    
    // Unit conversions 
    const double t_per_bin = 1.0;                       // Time (in above units) per bin, e.g., 1 ms per bin
    const double sample_rate = 1e4;                     // Sample rate (in above units), e.g., 1000 Hz
    
    // Data fields
    MatrixXd trial_data;                                // NxM matrix of doubles, rows as recording times (in "unit_time"), columns as trials, data values in "unit_data"
    MatrixXd spike_raster;                              // Nx2 matrix, each row one spike, columns as time (in "unit_time") and trial number
    
    // Analysis fields
    VectorXd autocorr;                                  // Autocorrelation of trial_data

  public:
    
    // Constructor and Destructor
    neuron(
      const int id_num = 0, 
      const std::string recording_name = "not_provided", 
      const std::string type = "generic", 
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
    
    // Member functions for fetching data
    MatrixXd fetch_trial_data() const;
    NumericMatrix fetch_trial_data_R() const;
    MatrixXd fetch_spike_raster() const;
    NumericMatrix fetch_spike_raster_R() const;
    
    // Member functions for fetching analysis results
    VectorXd fetch_autocorr() const;
    NumericVector fetch_autocorr_R() const;
    
    // Member functions for data analysis
    void compute_autocorrelation();

};

#endif
