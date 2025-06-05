#!/bin/bash
#SBATCH --job-name=parallel_test_mpirun
#SBATCH --output=parallel_mpirun_%j.out
#SBATCH --error=parallel_mpirun_%j.err
#SBATCH --nodes=2
#SBATCH --ntasks=80
#SBATCH --ntasks-per-node=40
#SBATCH --time=0-01:00:00

# Purge modules and load the required environment.
module purge
module load StdEnv/2023 gcc/13.3 openmpi/5.0.3

# Loop over the desired process counts.
for P in 1 8 20 40 80; do
  echo "======================================="
  echo "Running parallel read with P=$P using mpirun"
  echo "---------------------------------------"
  time mpirun -np $P ./parallel 1.1 morestepnumbers_1.7GB.dat 100000
  echo ""
done
