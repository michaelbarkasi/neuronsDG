# Fetch a neuron's identifying metadata

Fetch a neuron's identifying metadata

## Usage

``` r
fetch.id.data(nrn)
```

## Arguments

- nrn:

  A neuron object.

## Value

A list of the neuron's identifying fields and metadata (id number,
recording name, type, genotype, sex, hemisphere, region, age, whether
simulated, units, bin size, and sample rate).

## Examples

``` r
nrn <- new.neuron(id_num = 1, recording_name = "demo_recording", genotype = "WT")
fetch.id.data(nrn)
#> $id_num
#> [1] 1
#> 
#> $recording_name
#> [1] "demo_recording"
#> 
#> $type
#> [1] "generic"
#> 
#> $genotype
#> [1] "WT"
#> 
#> $sex
#> [1] "not_provided"
#> 
#> $hemi
#> [1] "not_provided"
#> 
#> $region
#> [1] "not_provided"
#> 
#> $age
#> [1] "not_provided"
#> 
#> $sim
#> [1] FALSE
#> 
#> $unit_time
#> [1] "ms"
#> 
#> $unit_sample_rate
#> [1] "Hz"
#> 
#> $unit_data
#> [1] "mV"
#> 
#> $t_per_bin
#> [1] 10
#> 
#> $sample_rate
#> [1] 10000
#> 
```
