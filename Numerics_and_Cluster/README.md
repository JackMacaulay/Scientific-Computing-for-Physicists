# PHY1610 Assignment 4: Numerics and Teach Cluster

## Overview
Simulating a marble's trajectory in a viscous fluid and determining friction rate.

## Build Instructions
- Generate test model data: `make testmodel.dat`
- Data analysis: `make analysis`

## Documentation
Generate documentation using Doxygen by running: make doc

## Job Script
SLURM job script (run_analysis.sh) is provided to run the analysis. Submit the job with: sbatch run_analysis.sh

## Robustness
Current frictionrate can give nan or inf values when delta v is very small. Some ideas for improvement:

Using epsilon threshold to avoid division by near-zero values
Weighted averaging to favour more reliable data
Adding error handling 