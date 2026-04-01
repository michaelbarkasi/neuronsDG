# Fit exponential decay function to autocorrelation of single neuron object

This function fits an exponential decay function to the autocorrelation
of a single neuron object using specified parameters. It's a wrapper
around the `fit_autocorrelation` method of the neuron class.

## Usage

``` r
fit.edf.autocorr(nrn, A0 = 0.001, tau0 = 1, ctol = 1e-08, max_evals = 500)
```

## Arguments

- nrn:

  Neuron object for which to fit autocorrelation.

- A0:

  Initial guess for amplitude parameter of exponential decay function
  (default: 0.001).

- tau0:

  Initial guess for time constant parameter of exponential decay
  function (default: 1.0).

- ctol:

  Convergence tolerance for fitting exponential decay function (default:
  1e-8).

- max_evals:

  Maximum number of evaluations for fitting exponential decay function
  (default: 500).

## Value

None. Modifies neuron in place.
