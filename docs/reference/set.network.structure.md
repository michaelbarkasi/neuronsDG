# Set network structure

This function sets the structure of a network object, defining its
layers, columns, neuron types, and local connectivity parameters. It
also generates local nodes based on the specified structure.

## Usage

``` r
set.network.structure(
  network,
  neuron_types = c("principal"),
  layer_names = c("layer"),
  n_layers = 1,
  n_columns = 1,
  layer_height = 250,
  column_width = 130,
  layer_separation_factor = 3,
  column_separation_factor = 3.5,
  neurons_per_node = 30,
  recurrence_factors = 0.5,
  pruning_threshold_factor = 0.1
)
```

## Arguments

- network:

  Network object to configure.

- neuron_types:

  Character vector giving types of neurons in the network, e.g.
  c("principal", "interneuron").

- layer_names:

  Character vector giving names of layers in the network, e.g. c("L2/3",
  "L4", "L5", "L6").

- n_layers:

  Integer giving number of layers in the network.

- n_columns:

  Integer giving number of columns in the network.

- layer_height:

  Numeric giving height of each layer (in units specified at network
  creation, default unit is microns, default value is 250.0).

- column_width:

  Numeric giving width of each column (in units specified at network
  creation, default unit is microns, default value is 130.0).

- layer_separation_factor:

  Numeric giving mean distance between layers as a fraction of layer
  height (default: 3.0).

- column_separation_factor:

  Numeric giving mean distance between columns as a fraction of column
  width (default: 3.5).

- neurons_per_node:

  Matrix giving number of neurons of each type per node in each layer;
  dimensions must match n_layers (rows) and length of neuron_types
  (columns).

- recurrence_factors:

  List of matrices giving local recurrence factors for each layer; each
  matrix must have dimensions matching length of neuron_types (rows and
  columns).

- pruning_threshold_factor:

  Numeric giving factor for pruning weak connections within nodes;
  connections with strength below this factor times the maximum
  connection strength in the node will be pruned (default: 0.1).

- neuron_type_valences:

  Numeric vector giving valences of each neuron type, e.g. c(1, -1) for
  excitatory and inhibitory neurons.

- neuron_type_temporal_modulation:

  Numeric matrix giving temporal modulation time components (for
  modulation time in the unit_time of the network) for each neuron type:
  bias, step size, and count cutoff (rows as neuron types, columns as
  components). Will example a single value or a vector of length three.

## Value

The updated network object with the specified structure and local nodes
generated.
