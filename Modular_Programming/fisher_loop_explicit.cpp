/******************************************************
 * Part 1: KPP-Fisher Solver with Explicit Triple Loop
 *
 * Solves:    u_t - u_xx = u*(1 - u)
 * BCs:       u(0,t)=A*sin²(t),  u(L,t)=0
 * IC:        u(x,0)=0  for x in [0,L]
 *
 * We form a full (N-2)x(N-2) matrix (for the interior points)
 * and multiply it by the interior of u using a loop.
 * Forward Euler is used for time stepping.
 *
 * Compile with:
 *   g++ -O2 -std=c++11 fisher_loop_explicit.cpp -o fisher_loop_explicit -lm
 ******************************************************/

 #include <iostream>    // for I/O
 #include <cmath>       // for sin, etc.
 #include "rarray.hpp"  // for rvector and rmatrix
 #include <string>      // for stoi, stod
 
 // Structure to store simulation parameters.
 struct SimulationParams {
     int P;       // Number of snapshots
     double L;    // Domain length
     double A;    // Amplitude at left boundary
     int N;       // Total grid points (incl. boundaries)
     double T;    // Final time
     double dt;   // Time step size
 };
 
 // Parse command-line arguments.
 // Expected: P L A N T dt
 SimulationParams parse_arguments(int argc, char* argv[]) {
     if (argc != 7) {
         std::cerr << "Usage: " << argv[0] << " P L A N T dt\n";
         std::exit(1);
     }
     SimulationParams params;
     params.P  = std::stoi(argv[1]);
     params.L  = std::stod(argv[2]);
     params.A  = std::stod(argv[3]);
     params.N  = std::stoi(argv[4]);
     params.T  = std::stod(argv[5]);
     params.dt = std::stod(argv[6]);
     if (params.N < 2) {
         std::cerr << "Error: N must be >= 2\n";
         std::exit(1);
     }
     return params;
 }
 
 // Build the full (N-2)x(N-2) matrix for the second derivative.
 // Only the interior points are used.
 rmatrix<double> build_diffusion_matrix(const SimulationParams& params) {
     int interiorCount = params.N - 2;
     double dx = params.L / (params.N - 1);
     rmatrix<double> M(interiorCount, interiorCount);
     double coeff = 1.0 / (dx * dx);
     for (int i = 0; i < interiorCount; i++) {
         M[i][i] = -2.0 * coeff;          // main diagonal
         if (i > 0)
             M[i][i - 1] = coeff;         // left neighbor
         if (i < interiorCount - 1)
             M[i][i + 1] = coeff;         // right neighbor
     }
     return M;
 }
 
 // Multiply the full matrix M by the interior of u_old using loops.
 // The result is stored in diff (size: N-2).
 void matvec_triple_loop(const rmatrix<double>& M,
                         const rvector<double>& u_old,
                         rvector<double>& diff) {
     int interiorCount = M.extent(0);
     // Clear diff vector
     for (int i = 0; i < interiorCount; i++)
         diff[i] = 0.0;
     // Multiply: note u_old index shifted by +1 for interior points.
     for (int i = 0; i < interiorCount; i++) {
         for (int j = 0; j < interiorCount; j++) {
             diff[i] += M[i][j] * u_old[j + 1];
         }
     }
 }
 
 // Print a snapshot: each line shows time, x, and u at that x.
 void print_snapshot(double t, const rvector<double>& u, double dx) {
     int N = u.size();
     for (int i = 0; i < N; i++)
         std::cout << t << " " << i * dx << " " << u[i] << "\n";
 }
 
 int main(int argc, char* argv[]) {
     // Get parameters from the command line.
     SimulationParams params = parse_arguments(argc, argv);
     double dx = params.L / (params.N - 1);
 
     // Initialize solution vectors (u_old and u_new).
     rvector<double> u_old(params.N), u_new(params.N);
     for (int i = 0; i < params.N; i++)
         u_old[i] = 0.0;
     // Set initial boundaries: u(0,0)=0, u(L,0)=0.
     u_old[0] = 0.0;
     u_old[params.N - 1] = 0.0;
 
     // Build the diffusion matrix for interior points.
     rmatrix<double> M = build_diffusion_matrix(params);
     int interiorCount = params.N - 2;
     rvector<double> diff(interiorCount);
 
     // Setup snapshot timing.
     double t = 0.0;
     int outputCount = 0;
     double outputInterval = (params.P > 1) ? params.T / (params.P - 1.0) : params.T;
     double nextOutputTime = 0.0;
 
     // Print initial snapshot.
     if (params.P > 0) {
         print_snapshot(t, u_old, dx);
         outputCount++;
         nextOutputTime += outputInterval;
     }
 
     // Main time-stepping loop (Forward Euler).
     while (t < params.T - 1e-14) {
         if (t + params.dt > params.T)
             params.dt = params.T - t;
 
         // Set boundary conditions at current time.
         u_old[0] = params.A * std::sin(t) * std::sin(t);
         u_old[params.N - 1] = 0.0;
 
         // Compute diffusion term for interior using our loop.
         matvec_triple_loop(M, u_old, diff);
 
         // Add contributions from boundaries to first and last interior nodes.
         if (interiorCount > 0) {
             double coeff = 1.0 / (dx * dx);
             diff[0] += coeff * u_old[0];
             diff[interiorCount - 1] += coeff * u_old[params.N - 1];
         }
 
         // Update interior points with Forward Euler (reaction term: u*(1-u)).
         for (int i = 0; i < interiorCount; i++) {
             int idx = i + 1; // shift index to match grid
             double u_val = u_old[idx];
             double reaction = u_val * (1.0 - u_val);
             u_new[idx] = u_val + params.dt * (diff[i] + reaction);
         }
 
         t += params.dt;
         // Enforce new boundary conditions.
         u_new[0] = params.A * std::sin(t) * std::sin(t);
         u_new[params.N - 1] = 0.0;
 
         u_old = u_new; // update for next time step
 
         // Print snapshot if it's time.
         if (outputCount < params.P && t >= nextOutputTime - 1e-14) {
             print_snapshot(t, u_old, dx);
             outputCount++;
             nextOutputTime += outputInterval;
         }
     }
     return 0;
 }
 