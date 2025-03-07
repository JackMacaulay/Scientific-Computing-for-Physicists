/******************************************************
 * fisher_sparse.cpp
 *
 * Part 3: KPP-Fisher Equation Solver with a Sparse
 *         Diffusion Matrix.
 *
 * PDE:   u_t - u_xx = u*(1 - u)
 * BCs:   u(0,t) = A*sin^2(t),  u(L,t) = 0
 * ICs:    u(x,0) = 0   for x in [0, L]
 *
 * The interior diffusion matrix is tridiagonal.
 * We only store nonzero entries to save time and memory.
 *
 ******************************************************/

 #include <iostream>
 #include <cmath>
 #include "rarray.hpp"
 #include <string>
 #include <vector>
 #include <cstdlib>
 
 // Holds simulation parameters.
 struct SimulationParams {
     int P;      // snapshots to print
     double L;   // domain length
     double A;   // amplitude at left boundary
     int N;      // grid points (incl. boundaries)
     double T;   // final time
     double dt;  // time step
 };
 
 // Stores one nonzero entry of the tridiagonal matrix.
 struct MatrixEntry {
     int row, col;
     double val;
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
 
 // Build the sparse tridiagonal matrix for the second derivative.
 std::vector<MatrixEntry> build_sparse_diffusion(const SimulationParams& p) {
     std::vector<MatrixEntry> spmat;
     int interiorCount = p.N - 2;
     spmat.reserve(interiorCount * 3);
     double dx = p.L / (p.N - 1);
     double coeff = 1.0 / (dx * dx);
     for (int i = 0; i < interiorCount; i++) {
         spmat.push_back({i, i, -2.0 * coeff});   // main diagonal
         if (i > 0) spmat.push_back({i, i - 1, coeff});  // left neighbor
         if (i < interiorCount - 1) spmat.push_back({i, i + 1, coeff}); // right neighbor
     }
     return spmat;
 }
 
 // Multiply the sparse matrix by the interior of u_old.
 void matvec_sparse(const std::vector<MatrixEntry>& spmat,
                    const rvector<double>& u_old,
                    rvector<double>& diff) {
     int interiorCount = diff.size();
     for (int i = 0; i < interiorCount; i++)
         diff[i] = 0.0;
     for (const auto& e : spmat)
         diff[e.row] += e.val * u_old[e.col + 1]; // shift for interior index
 }
 
 // Print one snapshot (time, x, u for every grid point).
 void print_snapshot(double t, const rvector<double>& u, double dx) {
     int N = u.size();
     for (int i = 0; i < N; i++)
         std::cout << t << " " << i * dx << " " << u[i] << "\n";
 }
 
 int main(int argc, char* argv[]) {
     // Read parameters.
     SimulationParams p = parse_arguments(argc, argv);
     double dx = p.L / (p.N - 1);
     
     // Set up solution vectors.
     rvector<double> u_old(p.N), u_new(p.N);
     for (int i = 0; i < p.N; i++) u_old[i] = 0.0;
     u_old[0] = 0.0; u_old[p.N - 1] = 0.0;
     
     // Build the sparse diffusion matrix.
     std::vector<MatrixEntry> spmat = build_sparse_diffusion(p);
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
     
     // Time stepping loop.
     while (t < p.T - 1e-14) {
         if (t + p.dt > p.T) p.dt = p.T - t;
         
         // Set boundary conditions.
         u_old[0] = p.A * std::sin(t) * std::sin(t);
         u_old[p.N - 1] = 0.0;
         
         // Diffusion term via sparse matrix multiply.
         matvec_sparse(spmat, u_old, diff);
         if (interiorCount > 0) {
             double coeff = 1.0 / (dx * dx);
             diff[0] += coeff * u_old[0];
             diff[interiorCount - 1] += coeff * u_old[p.N - 1];
         }
         
         // Update interior points with Forward Euler.
         for (int i_int = 0; i_int < interiorCount; i_int++) {
             int i_node = i_int + 1;
             double u_val = u_old[i_node];
             double reaction = u_val * (1.0 - u_val);
             u_new[i_node] = u_val + p.dt * (diff[i_int] + reaction);
         }
         
         t += p.dt;
         u_new[0] = p.A * std::sin(t) * std::sin(t);
         u_new[p.N - 1] = 0.0;
         u_old = u_new;
         
         if (outputCount < p.P && t >= nextOutputTime - 1e-14) {
             print_snapshot(t, u_old, dx);
             outputCount++;
             nextOutputTime += outputInterval;
         }
     }
     return 0;
 }
 