# Percolation Simulation (macaulay_a7)

This project simulates percolation in a 2D grid using random walkers.

## What It Does

### Lattice Setup:
- Creates a grid with **N** rows (height) and **M** columns (width).
- Each cell is "empty" with probability **p** or "filled" otherwise.

### Walker Simulation:
- Walkers start in each empty cell of the top row.
- They move left, right, up, or down:
  - Left/Right moves** have weight 1.0
  - Up moves have weight 1.0 / g
  - Down moves have weight g 
- Walkers move until they either reach the bottom row or exceed a maximum number of steps defined by:
  
  maxSteps = S * (N^2 + M^2)
  

### Results:
- The simulation counts how many walkers reach the bottom row and calculates the fraction of successful walkers.

## Files

- macaulay_a7.cpp
  Contains the complete simulation code.
- Makefile  
  Compiles the code and provides a run target to execute the simulation with default parameters.

## Requirements

- A C++ compiler supporting C++17 (e.g., g++).
- The rarray library 

## How to Build and Run

1. Build the executable:
   
   make
   
2. Run the simulation (with parameters: M=200, N=200, p=0.7, g=2.0, K=25, S=20):
   
   make run
   
3. Clean up the build files:
   
   make clean

## Notes

- A fixed random seed (`12345`) is used in the code for reproducibility.
- You can adjust the parameters to explore different percolation behaviors.

