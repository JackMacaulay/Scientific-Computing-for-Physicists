/******************************************************
 * Part 2: KPP-Fisher Solver with BLAS
 *
 * PDE:   u_t - u_xx = u*(1 - u)
 * BCs:   u(0,t) = A*sin²(t),  u(L,t) = 0
 * IC:    u(x,0) = 0 for x in [0, L]
 *
 * We use a full (N-2)x(N-2) matrix and call cblas_dgemv for
 * the diffusion term. Time stepping is done via Forward Euler.
 *
 ******************************************************/

 #include <iostream>
 #include <cmath>
 #include "rarray.hpp"
 #include <string>
 #include "cblas.h"
 
 // Holds simulation parameters.
 struct SimulationParams {
     int P;       // Number of snapshots
     double L;    // Domain length
     double A;    // Amplitude at left boundary
     int N;       // Number of grid points (including boundaries)
     double T;    // Final time
     double dt;   // Time step
 };
 
 // Parse command-line arguments.
 SimulationParams parse_arguments(int argc, char* argv[]) {
     if (argc != 7) {
         std::cerr << "Usage: " << argv[0] << " P L A N T dt\n";
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
         std::cerr << "Error: N must be >= 2\n";
         std::exit(1);
     }
     return p;
 }
 
 // Build the full diffusion matrix (for interior points).
 rmatrix<double> build_diffusion_matrix(const SimulationParams& p) {
     int interiorCount = p.N - 2;
     double dx = p.L / (p.N - 1);
     rmatrix<double> M(interiorCount, interiorCount);
     double coeff = 1.0 / (dx * dx);
     for (int i = 0; i < interiorCount; i++) {
         M[i][i] = -2.0 * coeff;
         if (i > 0)
             M[i][i-1] = coeff;
         if (i < interiorCount - 1)
             M[i][i+1] = coeff;
     }
     return M;
 }
 
 // Print one snapshot (time, x, u) for all grid points.
 void print_snapshot(double t, const rvector<double>& u, double dx) {
     int N = u.size();
     for (int i = 0; i < N; i++)
         std::cout << t << " " << i * dx << " " << u[i] << "\n";
 }
 
 int main(int argc, char* argv[]) {
     // Get parameters.
     SimulationParams p = parse_arguments(argc, argv);
     double dx = p.L / (p.N - 1);
 
     // Initialize solution arrays.
     rvector<double> u_old(p.N), u_new(p.N);
     for (int i = 0; i < p.N; i++)
         u_old[i] = 0.0;
     u_old[0] = 0.0; u_old[p.N-1] = 0.0;
 
     // Build diffusion matrix.
     rmatrix<double> M = build_diffusion_matrix(p);
     int interiorCount = p.N - 2;
     rvector<double> diff(interiorCount);
 
     // Setup snapshot timing.
     double t = 0.0;
     int outputCount = 0;
     double outputInterval = (p.P > 1) ? p.T / (p.P - 1.0) : p.T;
     double nextOutputTime = 0.0;
 
     // Print initial snapshot.
     if (p.P > 0) {
         print_snapshot(t, u_old, dx);
         outputCount++;
         nextOutputTime += outputInterval;
     }
 
     // Time-stepping loop.
     while (t < p.T - 1e-14) {
         if (t + p.dt > p.T)
             p.dt = p.T - t;
 
         // Update boundary conditions.
         u_old[0] = p.A * std::sin(t) * std::sin(t);
         u_old[p.N-1] = 0.0;
 
         // Compute diffusion: diff = M * (u_old interior)
         if (interiorCount > 0) {
             cblas_dgemv(CblasRowMajor, CblasNoTrans,
                         interiorCount, interiorCount,
                         1.0, &M[0][0], interiorCount,
                         &u_old[1], 1,
                         0.0, &diff[0], 1);
         }
 
         // Add boundary contributions.
         if (interiorCount > 0) {
             double coeff = 1.0 / (dx * dx);
             diff[0] += coeff * u_old[0];
             diff[interiorCount-1] += coeff * u_old[p.N-1];
         }
 
         // Forward Euler update for interior points.
         for (int i_int = 0; i_int < interiorCount; i_int++) {
             int i_node = i_int + 1;
             double u_val = u_old[i_node];
             double reaction = u_val * (1.0 - u_val);
             u_new[i_node] = u_val + p.dt * (diff[i_int] + reaction);
         }
 
         // Advance time.
         t += p.dt;
         u_new[0] = p.A * std::sin(t) * std::sin(t);
         u_new[p.N-1] = 0.0;
         u_old = u_new;
 
         // Print snapshot if it's time.
         if (outputCount < p.P && t >= nextOutputTime - 1e-14) {
             print_snapshot(t, u_old, dx);
             outputCount++;
             nextOutputTime += outputInterval;
         }
     }
 
     return 0;
 }
 