#OVERVIEW
There are two MPI implementations for processing a large dataset. Both programs compute the distribution of the logarithm (base 1.1) of the time step counts and output a normalized histogram.

1. Scatter - root process reads data in batches of size Z, then scatters each batch to the worker processes for histogram computation
2. Parallel - each MPI process reads its share of lines directly from the input file, skipping lines not owened by its rank.

#FILE STRUCTURE

scatter.cpp - MPI implementation where rank 0 reads data in batches
parallel.cpp - Contains MPI implementation where each rank reads its own subset of lines

Makefile - builds both scatter and parallel executables with make. Clean with make clean

run_scatter.sh/run_parallel.sh - slurm scripts for running each executable with different process count (P = 1, 8, 20, 40, 80)

#REQUIREMENTS
MPI library
C++11
The makefile will automatically load the relevant modules

#BUILD INSTRUCTIONS
Compile with make
Remove artifacts/executables with make clean

#USAGE
Run slurm scripts with

sbatch run_scatter.sh or
sbatch run_parallel.sh

These scripts will send a job to the cluster for either the parallel or scatter processes and loop over the process count P = 1, 8, 20, 40, 80, time each run and redirect the output histogram into files name scatter_jobnumber or parallel_jobnumber


#OUTPUT
Each executable will print a 30-bin histogram with two columns:

1. Bin start (in log domain)
2. Fraction of total data points in that bin


