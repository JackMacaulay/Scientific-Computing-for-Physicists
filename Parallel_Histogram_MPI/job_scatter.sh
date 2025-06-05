#!/bin/bash
#SBATCH --job-name=scatter_test_mpirun
#SBATCH --output=scatter_mpirun_%j.out
#SBATCH --error=scatter_mpirun_%j.err
#SBATCH --nodes=2                   # Request 2 nodes.
#SBATCH --ntasks=80                 # Allocate 80 tasks (max needed for P=80).
#SBATCH --ntasks-per-node=40        # Use 40 tasks per node.
#SBATCH --time=0-01:00:00           # Set the maximum job time to 1 hour.

# Purge modules and load the required environment.
module purge
module load StdEnv/2023 gcc/13.3 openmpi/5.0.3

# Loop over the desired process counts.
for P in 1 8 20 40 80; do
  echo "======================================="
  echo "Running scatter with P=$P using mpirun"
  echo "---------------------------------------"
  time mpirun -np $P ./scatter 1.1 morestepnumbers_1.7GB.dat 100000
  echo ""
done
