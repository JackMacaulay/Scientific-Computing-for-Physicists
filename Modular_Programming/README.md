# KPP-Fisher Equation Simulation

This project contains three C++ codes that solve the KPP-Fisher equation:

  uₜ - uₓₓ = u (1 - u)

with boundary conditions:
  u(0,t) = A · sin²(t)  and  u(L,t) = 0

and initial condition:
  u(x,0) = 0 for x ∈ [0, L].

There are three implementations:
- **fisher_loop_explicit.cpp** – Uses a full (N-2)x(N-2) matrix with explicit loops.
- **fisher_blas.cpp** – Uses a full matrix with BLAS (cblas_dgemv) for multiplication.
- **fisher_sparse.cpp** – Uses a sparse (tridiagonal) matrix representation.

## Building

Have OpenBLAS installed. Then, in the project directory, run:

```bash
make

## Running

Each executable requires six parameters:

P - Number of snapshots
L - Domain length
A - Amplitude at x=0
N - Number of grid points (boundary included)
T - final simulation time
dt = time step

e.g. formatted P L A N T
./fisher_loop_explicit 400 5.0 0.2 100 10.0 0.001
./fisher_blas         400 5.0 0.2 100 10.0 0.001
./fisher_sparse       400 5.0 0.2 100 10.0 0.001