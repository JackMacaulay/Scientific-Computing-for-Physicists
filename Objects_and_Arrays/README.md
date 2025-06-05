# KPP-Fisher Equation Solver (Assignment 6)

This repository contains a solver for the KPP-Fisher equation using an operator splitting method. The equation is

\[
\frac{\partial u}{\partial t} - \frac{\partial^2 u}{\partial x^2} = u(1 - u)
\]

with **periodic boundary conditions** \( u(0,t) = u(L,t) \) and the **initial condition**

\[
u(x,0) = A \left[\sin\left(\frac{\pi x}{L}\right)\right]^{100}.
\]

## Method

The solution is computed by splitting each time step into two substeps:

1. **Reaction Substep:**  
   Solves  
   \[
   \frac{\partial u}{\partial t} = u(1-u)
   \]
   exactly over a time interval \(\Delta t\) using the formula:
   \[
   u(x,t+\Delta t)=\frac{u(x,t)}{u(x,t)+[1-u(x,t)]e^{-\Delta t}}.
   \]
   
2. **Diffusion Substep:**  
   Solves  
   \[
   \frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}
   \]
   by transforming \(u\) to the Fourier domain using FFTW, updating each Fourier mode with:
   \[
   \hat{u}(k,t+\Delta t)=\hat{u}(k,t)e^{-k^2\Delta t},
   \]
   and then applying the inverse Fourier transform. Here, the wave numbers are given by \( k = \frac{2\pi q}{L} \).

The operator splitting is applied repeatedly for a total simulation time \(T\).

## Parameters

The simulation uses the following parameters (passed as command-line arguments):

- **P**: Number of snapshots to output
- **L**: Domain length (set to 5)
- **A**: Amplitude in the initial condition (set to 1)
- **N**: Number of grid points (set to 100)
- **T**: Final simulation time (set to 10)
- **dt**: Time step (set to 0.01)

## Files

- `macaulay_a6.cpp`: The C++ source code for the solver.
- `Makefile`: A makefile to compile the source code.
- `README.md`: This readme file.

## Compilation

Ensure that the `gcc`, `fftw`, and `rarray` modules are loaded. Then compile using:
```bash
make
