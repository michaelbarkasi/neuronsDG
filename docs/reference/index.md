# Package index

## Initialization functions

- [`new.motif()`](https://michaelbarkasi.github.io/neurons/reference/new.motif.md)
  : Initialize network (circuit) motif
- [`new.network()`](https://michaelbarkasi.github.io/neurons/reference/new.network.md)
  : Initialize neuron network
- [`new.neuron()`](https://michaelbarkasi.github.io/neurons/reference/new.neuron.md)
  : Initialize neuron
- [`load.rasters.as.neurons()`](https://michaelbarkasi.github.io/neurons/reference/load.rasters.as.neurons.md)
  : Load neuron data from raster file

## Neuron method wrappers

- [`fit.edf.autocorr()`](https://michaelbarkasi.github.io/neurons/reference/fit.edf.autocorr.md)
  : Fit exponential decay function to autocorrelation of single neuron
  object
- [`compute.autocorr()`](https://michaelbarkasi.github.io/neurons/reference/compute.autocorr.md)
  : Compute autocorrelation of single neuron object
- [`estimate.autocorr.params()`](https://michaelbarkasi.github.io/neurons/reference/estimate.autocorr.params.md)
  : Function to estimate autocorrelation parameters using dichotomized
  Gaussian simulations
- [`process.autocorr()`](https://michaelbarkasi.github.io/neurons/reference/process.autocorr.md)
  : Function to process autocorrelation of neuron list
- [`plot.raster()`](https://michaelbarkasi.github.io/neurons/reference/plot-raster.md)
  : Function to plot spike raster of neuron
- [`plot.autocorrelation()`](https://michaelbarkasi.github.io/neurons/reference/plot-autocorrelation.md)
  : Function to plot autocorrelation of neuron

## Motif method wrappers

- [`load.projection.into.motif()`](https://michaelbarkasi.github.io/neurons/reference/load.projection.into.motif.md)
  : Load projection into motif

## Network method wrappers

- [`set.network.structure()`](https://michaelbarkasi.github.io/neurons/reference/set.network.structure.md)
  : Set network structure
- [`apply.circuit.motif()`](https://michaelbarkasi.github.io/neurons/reference/apply.circuit.motif.md)
  : Apply circuit motif to network
- [`run.SGT()`](https://michaelbarkasi.github.io/neurons/reference/run.SGT.md)
  : Run Spatial Growth-Transform network simulation
- [`plot.network()`](https://michaelbarkasi.github.io/neurons/reference/plot-network.md)
  : Plot network as directed graph
- [`plot.network.traces()`](https://michaelbarkasi.github.io/neurons/reference/plot-network-traces.md)
  : Plot spike traces for network from SGT simulation

## Network cell type functions

- [`print.known.celltypes()`](https://michaelbarkasi.github.io/neurons/reference/print-known-celltypes.md)
  : Print known cell types
- [`fetch.cell.type.params()`](https://michaelbarkasi.github.io/neurons/reference/fetch.cell.type.params.md)
  : Fetch cell type parameters
- [`add.cell.type()`](https://michaelbarkasi.github.io/neurons/reference/add.cell.type.md)
  : Add new cell type
- [`modify.cell.type()`](https://michaelbarkasi.github.io/neurons/reference/modify.cell.type.md)
  : Modify existing cell type

## Analysis functions

- [`analyze.autocorr()`](https://michaelbarkasi.github.io/neurons/reference/analyze.autocorr.md)
  : Function to analyze autocorrelation estimates with bootstrapping
- [`test.sigma.assumption()`](https://michaelbarkasi.github.io/neurons/reference/test.sigma.assumption.md)
  : Helper function for testing assumption of autocorrelation processing

## Preprocessing functions

- [`import.kilo4()`](https://michaelbarkasi.github.io/neurons/reference/import.kilo4.md)
  : Import raw kilosort4 data
- [`preprocess.kilo4()`](https://michaelbarkasi.github.io/neurons/reference/preprocess.kilo4.md)
  : Preprocess kilosort4 data
- [`summarize.cluster.key()`](https://michaelbarkasi.github.io/neurons/reference/summarize.cluster.key.md)
  : Summarize cluster key by covariates
