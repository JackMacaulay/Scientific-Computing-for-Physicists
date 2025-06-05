SHELL := /bin/bash

MPICXX = mpicxx
CXXFLAGS = -O3 -std=c++11

# Build both versions by default.
all: scatter parallel

# "scatter" version (root reads batches & scatters)
scatter: scatter.cpp
	/bin/bash -l -c "module load StdEnv/2023 gcc/13.3 openmpi/5.0.3 && \
	    $(MPICXX) $(CXXFLAGS) scatter.cpp -o scatter"

# "parallel" version (step 5: each rank reads subset in parallel)
parallel: parallel.cpp
	/bin/bash -l -c "module load StdEnv/2023 gcc/13.3 openmpi/5.0.3 && \
	    $(MPICXX) $(CXXFLAGS) parallel.cpp -o parallel"

clean:
	rm -f scatter parallel

