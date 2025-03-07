/******************************************************
 * fisher_sparse.cpp
 *
 * Part 3: KPP-Fisher Equation Solver using a Sparse
 *         Matrix Representation for the Diffusion Operator.
 *
 * The PDE is:
 *     u_t - u_xx = u*(1 - u)
 *
 * with boundary conditions:
 *     u(0,t) = A*sin^2(t),    u(L,t) = 0,
 *
 * and initial condition:
 *     u(x,0) = 0  for x in [0,L].
 *
 * The domain is discretized into N points, and the
 * solution is evolved in time using Forward Euler
 * with time step Δt (chosen for numerical stability).
 *
 * Only the nonzero elements of the (N-2)x(N-2) interior
 * diffusion matrix (which is tridiagonal) are stored.
 * This saves memory and computation, since there are only
 * ~3 nonzeros per row.
 *
 * The simulation prints P snapshots. Each snapshot prints
 * all grid points (with each line containing: time, x, u(x,t)).
 *
 * Compile with:
 *     g++ -O2 -std=c++11 fisher_sparse.cpp -o fisher_sparse -lm
 ******************************************************/

#include <iostream>    // for std::cout, std::cerr
#include <cmath>       // for sin, etc.
#include "rarray.hpp"      // for rvector (1D) container
#include <string>      // for stoi, stod
#include <vector>      // for std::vector to hold sparse entries
#include <cstdlib>     // for std::exit

// Structure to hold simulation parameters
struct SimulationParams {
    int    P;    // number of snapshots to output
    double L;    // domain length
    double A;    // amplitude for left boundary (u(0,t))
    int    N;    // number of grid points (including boundaries)
    double T;    // final time
    double dt;   // time step size
};

// Structure to hold a nonzero entry of the diffusion matrix
struct MatrixEntry {
    int row;    // interior row index (0 <= row < N-2)
    int col;    // interior column index (0 <= col < N-2)
    double val; // value of the matrix element
};

//-----------------------------------------------------------
// Function: parse_arguments
// Purpose:  Parse command-line arguments into a SimulationParams struct.
// Expected arguments: P L A N T dt
//-----------------------------------------------------------
SimulationParams parse_arguments(int argc, char* argv[]) {
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " P L A N T dt\n"
                  << "  P   = number of snapshots\n"
                  << "  L   = domain length\n"
                  << "  A   = amplitude for u(0,t)\n"
                  << "  N   = number of grid points (including boundaries)\n"
                  << "  T   = final time\n"
                  << "  dt  = time step size\n"
                  << "Example:\n  " << argv[0] << " 400 5.0 0.2 100 10 0.001\n";
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
        std::cerr << "Error: N must be >= 2 (including boundaries)\n";
        std::exit(1);
    }
    return p;
}

//-----------------------------------------------------------
// Function: build_sparse_diffusion
// Purpose:  Build a sparse representation of the (N-2)x(N-2)
//           diffusion matrix corresponding to the second
//           derivative (finite difference) operator.
//           For each interior point i (0-based index for interior),
//           the nonzero entries are:
//             - Main diagonal: (i, i) = -2/(dx^2)
//             - Sub-diagonal:   (i, i-1) =  1/(dx^2)   if i > 0
//             - Super-diagonal: (i, i+1) =  1/(dx^2)   if i < interiorCount-1
//-----------------------------------------------------------
std::vector<MatrixEntry> build_sparse_diffusion(const SimulationParams& p) {
    std::vector<MatrixEntry> spmat;
    int interiorCount = p.N - 2;
    spmat.reserve(interiorCount * 3); // reserve space for about 3 nonzeros per row

    double dx = p.L / (p.N - 1);
    double coeff = 1.0 / (dx * dx);

    // Loop over each interior row (i_int corresponds to actual index i = i_int+1)
    for (int i_int = 0; i_int < interiorCount; i_int++) {
        // Main diagonal entry
        spmat.push_back({i_int, i_int, -2.0 * coeff});
        // Sub-diagonal entry (if exists)
        if (i_int > 0) {
            spmat.push_back({i_int, i_int - 1, coeff});
        }
        // Super-diagonal entry (if exists)
        if (i_int < interiorCount - 1) {
            spmat.push_back({i_int, i_int + 1, coeff});
        }
    }
    return spmat;
}

//-----------------------------------------------------------
// Function: matvec_sparse
// Purpose:  Compute the product of the sparse matrix (given by spmat)
//           and the interior portion of the solution vector u_old.
//           The interior portion is u_old[1] ... u_old[N-2] (since boundaries
//           are at indices 0 and N-1).
//           The result is stored in diff (of size N-2).
//-----------------------------------------------------------
void matvec_sparse(const std::vector<MatrixEntry>& spmat,
                   const rvector<double>& u_old,
                   rvector<double>& diff) {
    int interiorCount = diff.size();
    // Zero out the result vector
    for (int i = 0; i < interiorCount; i++) {
        diff[i] = 0.0;
    }
    // For each nonzero entry, accumulate its contribution.
    // Note: Since spmat uses interior indices, we use u_old[entry.col + 1] to access the corresponding grid point.
    for (const auto& entry : spmat) {
        diff[entry.row] += entry.val * u_old[entry.col + 1];
    }
}

//-----------------------------------------------------------
// Function: print_snapshot
// Purpose:  Print the current snapshot of the solution.
//           Each line prints: time, x-position, u-value for every grid point.
//-----------------------------------------------------------
void print_snapshot(double t, const rvector<double>& u, double dx) {
    int N = u.size();
    for (int i = 0; i < N; i++) {
        double x = i * dx;
        std::cout << t << " " << x << " " << u[i] << "\n";
    }
}

//-----------------------------------------------------------
// Main: Solve the KPP-Fisher Equation using a Sparse Matrix Approach
//-----------------------------------------------------------
int main(int argc, char* argv[]) {
    // Parse input parameters from the command line.
    SimulationParams p = parse_arguments(argc, argv);
    
    // Compute spatial step size.
    double dx = p.L / (p.N - 1);

    // Allocate solution vectors (u_old for current time, u_new for next time).
    rvector<double> u_old(p.N), u_new(p.N);
    // Initialize u(x,0) = 0 for all grid points.
    for (int i = 0; i < p.N; i++) {
        u_old[i] = 0.0;
    }
    // Set the boundary conditions at t = 0.
    u_old[0] = 0.0;         // u(0,0) = A*sin^2(0) = 0
    u_old[p.N - 1] = 0.0;     // u(L,0) = 0

    // Build the sparse representation of the diffusion matrix.
    std::vector<MatrixEntry> spmat = build_sparse_diffusion(p);
    int interiorCount = p.N - 2; // number of interior grid points

    // Allocate vector for the result of the sparse matrix-vector product.
    rvector<double> diff(interiorCount);

    // Set up snapshot scheduling.
    double t = 0.0;                    // current time
    int outputCount = 0;               // count of snapshots printed
    double outputInterval = (p.P > 1) ? p.T / (p.P - 1.0) : p.T;
    double nextOutputTime = 0.0;

    // Print the initial snapshot (t = 0).
    if (p.P > 0) {
        print_snapshot(t, u_old, dx);
        outputCount++;
        nextOutputTime += outputInterval;
    }

    // Main time-stepping loop (Forward Euler).
    while (t < p.T - 1e-14) {
        // Adjust dt if necessary so as not to overshoot final time.
        if (t + p.dt > p.T) {
            p.dt = p.T - t;
        }

        // Update boundary conditions at the current time.
        u_old[0] = p.A * std::sin(t) * std::sin(t);  // Left boundary: u(0,t) = A*sin^2(t)
        u_old[p.N - 1] = 0.0;                         // Right boundary: u(L,t) = 0

        // Compute the diffusion term for interior points:
        // diff = M * (u_old interior)
        // (using our sparse representation)
        matvec_sparse(spmat, u_old, diff);

        // Add contributions from boundary neighbors:
        // The first interior point (u_old[1]) has a left neighbor u_old[0],
        // and the last interior point (u_old[N-2]) has a right neighbor u_old[N-1].
        if (interiorCount > 0) {
            double coeff = 1.0 / (dx * dx);
            diff[0] += coeff * u_old[0];
            diff[interiorCount - 1] += coeff * u_old[p.N - 1];
        }

        // Update interior points with the Forward Euler method:
        // u_new[i] = u_old[i] + dt*(diffusion term + reaction term)
        for (int i_int = 0; i_int < interiorCount; i_int++) {
            int i_node = i_int + 1;  // shift to actual grid index
            double u_val = u_old[i_node];
            double reaction = u_val * (1.0 - u_val); // reaction term: u(1 - u)
            u_new[i_node] = u_val + p.dt * (diff[i_int] + reaction);
        }

        // Advance time.
        t += p.dt;

        // Enforce boundary conditions at new time.
        u_new[0] = p.A * std::sin(t) * std::sin(t);
        u_new[p.N - 1] = 0.0;

        // Update u_old with the new solution.
        u_old = u_new;

        // If it is time for a snapshot, print the current solution.
        if (outputCount < p.P && t >= nextOutputTime - 1e-14) {
            print_snapshot(t, u_old, dx);
            outputCount++;
            nextOutputTime += outputInterval;
        }
    }

    return 0;
}
