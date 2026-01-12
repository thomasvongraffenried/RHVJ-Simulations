# RHV-J Simulations
This repository contains four custom R scripts used to simulate and post-process the spatial evolution of a non-recombining viral population in a continuous two-dimensional landscape. The simulations were developed to investigate whether occasional long-distance vector-mediated transmission events can generate phylogenetic patterns similar to those observed in empirical RHV-J data. 

The two scripts differ only in the dispersal kernel used for the placement of the offspring. 

## Repository Contents

- LDT_Seeds1to27.R

  Simulates viral evolution with both short-range and rare long-distance transmission events. The simulation is repeated for 27 independent random seeds. 

- NoLongDistanceDispersal_Simulation.R

  Simulates viral evolution with only short-range transmission events. The simulation is repeated for 27 independent random seeds.

- LDT_FindSuitableGeneration.R

  Post-processes the simulation outputs for each of the 27 independent seeds to identify the earliest generation with sufficient individuals in predefined spatial areas, and performs spatially balanced sampling of sequences.

- nLDT_FindSuitableGeneration.R
 
  Post-processes the simulation outputs for each of the 27 independent seeds to identify the earliest generation with sufficient individuals in predefined spatial areas, and performs spatially balanced sampling of sequences.

- README.md

  Description of simulation framework and usage.

## Outputs
The simulation and post-processing scripts generate:

- FASTA files containing viral sequences for each seed and generation
- Text files with spatial coordinates of individual samples
- Sampled FASTA and NEXUS files for predefined spatial areas
- Optional spatial plots illustrating population distributions and sampling areas

## Requirements 
R (≥ 4.0) and standard CRAN packages listed at the top of each script.
