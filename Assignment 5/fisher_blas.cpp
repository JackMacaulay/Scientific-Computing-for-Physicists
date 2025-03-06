/******************************************************
 * Part 2: KPP-Fisher equation solved with a full (N-2)x(N-2)
 * matrix but using a BLAS routine (cblas_dgemv) for
 * the matrix-vector multiplication. Forward Euler in time.
 *
 * Compile with:
 *   g++ -O2 -std=c++11 fisher_blas.cpp -o fisher_blas -lopenblas -lm
 ******************************************************/

#include <iostream>    // for std::cout, std::cerr
#include <cmath>       // for sin, etc.
#include "rarray.hpp"      // for rvector, rmatrix
#include <string>      // for stoi, stod
#include "cblas.h"

// A struct to hold all simulation parameters
struct SimulationParams {
    int P;        // number of snapshots
    double L;     // domain length
    double A;     // amplitude for left boundary
    int N;        // number of grid points (including boundaries)
    double T;     // final time
    double dt;    // time step
};

//------------------------------------------------------------------------------
// 1) Function to parse the command-line arguments
//    We expect: P L A N T dt
//------------------------------------------------------------------------------
SimulationParams parse_arguments(int argc, char* argv[])
{
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0]
                  << " P L A N T dt\n\n"
                  << "  P   = number of snapshots\n"
                  << "  L   = domain length\n"
                  << "  A   = amplitude at x=0 boundary\n"
                  << "  N   = number of grid points\n"
                  << "  T   = final time\n"
                  << "  dt  = time-step size\n\n"
                  << "Example:\n"
                  << "  " << argv[0] << " 400 5.0 0.2 100 10.0 0.001\n"
                  << std::endl;
        std::exit(1);
    }

    SimulationParams p;
    p.P  = std::stoi(argv[1]);
    p.L  = std::stod(argv[2]);
    p.A  = std::stod(argv[3]);
    p.N  = std::stoi(argv[4]);
    p.T  = std::stod(argv[5]);
    p.dt = std::stod(argv[6]);

    if (p.N < 2) {
        std::cerr << "Error: N must be >= 2 (including boundary points)\n";
        std::exit(1);
    }
    return p;
}

//------------------------------------------------------------------------------
// 2) Function to build an (N-2)x(N-2) diffusion matrix for the second derivative
//    in 1D with spacing dx. We'll treat the boundary conditions separately.
//------------------------------------------------------------------------------
rmatrix<double> build_diffusion_matrix(const SimulationParams& p)
{
    int interiorCount = p.N - 2;
    double dx = p.L / (p.N - 1);

    rmatrix<double> M(interiorCount, interiorCount);

    // Initialize everything to 0
    for (int i = 0; i < interiorCount; i++) {
        for (int j = 0; j < interiorCount; j++) {
            M[i][j] = 0.0;
        }
    }

    double coeff = 1.0 / (dx * dx);
    for (int i = 0; i < interiorCount; i++) {
        // main diagonal
        M[i][i] = -2.0 * coeff;

        // sub-diagonal
        if (i > 0) {
            M[i][i - 1] = coeff;
        }
        // super-diagonal
        if (i < interiorCount - 1) {
            M[i][i + 1] = coeff;
        }
    }
    return M;
}

//------------------------------------------------------------------------------
// 3) Function to print a snapshot of the solution for all N points
//    Format: time, position, u-value
//------------------------------------------------------------------------------
void print_snapshot(double t, const rvector<double>& u, double dx)
{
    int N = u.size();
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        std::cout << t << " " << x << " " << u[i] << "\n";
    }
}

//------------------------------------------------------------------------------
// 4) Main program: solve the KPP-Fisher equation using BLAS for M*u
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // 4.1 parse parameters
    SimulationParams p = parse_arguments(argc, argv);

    // 4.2 allocate solution arrays (old and new)
    rvector<double> u_old(p.N), u_new(p.N);
    for (int i = 0; i < p.N; i++) {
        u_old[i] = 0.0;  // initial condition: u(x,0) = 0
    }
    // boundary at t=0: left BC = 0.0, right BC = 0.0
    u_old[0]       = 0.0;
    u_old[p.N-1]   = 0.0;

    // 4.3 build the (N-2)x(N-2) matrix for the diffusion
    rmatrix<double> M = build_diffusion_matrix(p);
    int interiorCount  = p.N - 2;

    // vector to store M*u_interior
    rvector<double> diff(interiorCount);

    double dx = p.L / (p.N - 1);

    // 4.4 scheduling for output
    double t = 0.0;
    int outputCount = 0;
    double outputInterval = (p.P > 1) ? p.T / (p.P - 1.0) : p.T;
    double nextOutputTime = 0.0;

    // print initial snapshot if P>0
    if (p.P > 0) {
        print_snapshot(t, u_old, dx);
        outputCount++;
        nextOutputTime += outputInterval;
    }

    // 4.5 main time-stepping loop
    while (t < p.T - 1e-14) {
        // don't overshoot final time
        if (t + p.dt > p.T) {
            p.dt = p.T - t;
        }

        // update boundary at current time
        u_old[0]             = p.A * std::sin(t) * std::sin(t);  // left boundary
        u_old[p.N - 1]       = 0.0;                              // right boundary

        // (a) use BLAS: diff = M*(u_interior)
        // cblas_dgemv: y = alpha*M*x + beta*y
        // here: alpha=1, beta=0, incX=1, incY=1
        if (interiorCount > 0) {
            cblas_dgemv(
                CblasRowMajor,      // row-major storage
                CblasNoTrans,       // no transpose
                interiorCount,      // # rows
                interiorCount,      // # cols
                1.0,                // alpha
                &M[0][0],           // pointer to matrix data
                interiorCount,      // leading dimension = # of cols for row-major
                &u_old[1],          // pointer to interior of u_old (skip boundary)
                1,                  // incX
                0.0,                // beta
                &diff[0],           // output
                1                   // incY
            );
        }

        // (b) Add boundary neighbor contributions
        if (interiorCount > 0) {
            double coeff = 1.0 / (dx*dx);
            diff[0]                += coeff * u_old[0];
            diff[interiorCount-1]  += coeff * u_old[p.N - 1];
        }

        // (c) forward euler update for interior
        for (int i_int = 0; i_int < interiorCount; i_int++) {
            int i_node = i_int + 1;  // shift by +1 to get the actual domain index
            double val_old  = u_old[i_node];
            double reaction = val_old * (1.0 - val_old);  // KPP-Fisher: u(1-u)
            u_new[i_node]   = val_old + p.dt * (diff[i_int] + reaction);
        }

        // advance time
        t += p.dt;

        // boundary at new time
        u_new[0]             = p.A * std::sin(t) * std::sin(t);
        u_new[p.N - 1]       = 0.0;

        // swap
        u_old = u_new;

        // output snapshots if time >= nextOutputTime
        if (outputCount < p.P && t >= nextOutputTime - 1e-14) {
            print_snapshot(t, u_old, dx);
            outputCount++;
            nextOutputTime += outputInterval;
        }
    }

    return 0;
}
