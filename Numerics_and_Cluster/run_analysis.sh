#!/bin/bash
#SBATCH --job-name=assignment4_analysis
#SBATCH --output=analysis_output.txt
#SBATCH --error=analysis_error.txt
#SBATCH --time=00:10:00
#SBATCH --partition=compute

module load StdEnv/2023
module load gcc/12.3
module load rarray/2.8.0
module load boost/1.85.0

make
./analyze -f testmodel.dat -o analysis.out
cat analysis.out
make clean
