#include <iostream>    // for std::cout, std::cerr
#include <cmath>       // for sin, etc.
#include <rarray>      // for rvector and rmatrix
#include <string>      // for stoi, stod

// A simple struct to store all relevant parameters.
struct SimulationParams {
    int P;            // Number of snapshots to output
    double L;         // Domain length
    double A;         // Amplitude for boundary condition at x=0
    int N;            // Number of grid points (including boundaries)
    double T;         // Final time
    double dt;        // Time step
};

//------------------------------------------------------------------------------
// 1) Parse command-line arguments and store in SimulationParams.
//    We expect: P L A N T dt
//------------------------------------------------------------------------------
SimulationParams parse_arguments(int argc, char* argv[])
{
    if (argc != 7) {
        std::cerr 
            << "Usage: " << argv[0] 
            << " P L A N T dt\n\n"
            << "  P   = number of snapshots\n"
            << "  L   = domain length\n"
            << "  A   = amplitude for BC at x=0\n"
            << "  N   = number of grid points (including boundaries)\n"
            << "  T   = final time\n"
            << "  dt  = time-step size\n\n"
            << "Example:\n"
            << "  " << argv[0] << " 400 5.0 0.2 100 10.0 0.001\n"
            << std::endl;
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
        std::cerr << "Error: N must be >= 2 (including boundary points)\n";
        std::exit(1);
    }
    return params;
}

//------------------------------------------------------------------------------
// 2) Build the (N-2)x(N-2) diffusion matrix for the 1D second-derivative operator.
//    We'll treat boundary conditions separately, so this matrix only applies to
//    the interior points.
//------------------------------------------------------------------------------
rmatrix<double> build_diffusion_matrix(const SimulationParams& params)
{
    int interiorCount = params.N - 2;
    double dx = params.L / (params.N - 1);

    rmatrix<double> M(interiorCount, interiorCount);

    // Initialize all entries to zero
    for (int i = 0; i < interiorCount; i++) {
        for (int j = 0; j < interiorCount; j++) {
            M[i][j] = 0.0;
        }
    }

    // Standard second-difference coefficients: -2 on diagonal, +1 on sub & super
    double coeff = 1.0 / (dx * dx);
    for (int i = 0; i < interiorCount; i++) {
        M[i][i] = -2.0 * coeff;  
        if (i > 0) {
            M[i][i - 1] = coeff; 
        }
        if (i < interiorCount - 1) {
            M[i][i + 1] = coeff;
        }
    }
    return M;
}

//------------------------------------------------------------------------------
// 3) Multiply the matrix M by the interior portion of u_old using an explicit
//    triple loop. The result is stored in 'diff', which has size (N-2).
//------------------------------------------------------------------------------
void matvec_triple_loop(const rmatrix<double>& M,
                        const rvector<double>& u_old,
                        rvector<double>& diff)
{
    int interiorCount = M.extent(0);

    // First, clear diff
    for (int i_int = 0; i_int < interiorCount; i_int++) {
        diff[i_int] = 0.0;
    }

    // Triple loop (row-by-row, col-by-col)
    // Actually it's a double loop, but it's “triple” if we consider
    // the overall code structure for time stepping as well.
    for (int i_int = 0; i_int < interiorCount; i_int++) {
        for (int j_int = 0; j_int < interiorCount; j_int++) {
            // shift by +1 to index the correct interior location in u_old
            diff[i_int] += M[i_int][j_int] * u_old[j_int + 1];
        }
    }
}

//------------------------------------------------------------------------------
// 4) Print the solution (time, position, value) for all grid points.
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
// 5) Main driver program
//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
    // Parse parameters from command line
    SimulationParams params = parse_arguments(argc, argv);

    // Compute dx based on N and L
    double dx = params.L / (params.N - 1);

    // Allocate vectors to store old and new solutions at each time step
    rvector<double> u_old(params.N), u_new(params.N);

    // Initial condition:  u(x,0) = 0
    for (int i = 0; i < params.N; i++) {
        u_old[i] = 0.0;
    }
    // Enforce boundary conditions at t=0
    // left boundary = A*sin^2(0) = 0, right boundary = 0
    u_old[0]            = 0.0;
    u_old[params.N-1]   = 0.0;

    // Build the (N-2)x(N-2) matrix
    rmatrix<double> M = build_diffusion_matrix(params);

    // We'll store M * (u_old interior) in this vector
    int interiorCount = params.N - 2;
    rvector<double> diff(interiorCount);

    // Time-stepping and output configuration
    double t = 0.0;             // current time
    int outputCount = 0;        // number of snapshots printed

    // Compute how often we print a snapshot
    double outputInterval;
    if (params.P > 1) {
        outputInterval = params.T / (params.P - 1.0);
    } else {
        outputInterval = params.T;
    }
    double nextOutputTime = 0.0;

    // Print initial snapshot if P>0
    if (params.P > 0) {
        print_snapshot(t, u_old, dx);
        outputCount++;
        nextOutputTime += outputInterval;
    }

    // Main time-stepping loop
    while (t < params.T - 1e-14) 
    {
        // If dt would overshoot T, adjust it
        if (t + params.dt > params.T) {
            params.dt = params.T - t;
        }

        // 1) Boundary condition at current time
        u_old[0]             = params.A * std::sin(t) * std::sin(t);
        u_old[params.N - 1]  = 0.0;

        // 2) Diffusion computation using the triple loop
        matvec_triple_loop(M, u_old, diff);

        // 3) Account for boundary neighbors for the first & last interior cell
        if (interiorCount > 0) {
            double coeff = 1.0 / (dx * dx);
            diff[0]               += coeff * u_old[0];
            diff[interiorCount-1] += coeff * u_old[params.N-1];
        }

        // 4) Update each interior node with forward Euler
        for (int i_int = 0; i_int < interiorCount; i_int++) {
            int i_node = i_int + 1;  // shift to actual domain index
            double val_old = u_old[i_node];
            double reaction = val_old * (1.0 - val_old); // KPP-Fisher: u(1-u)
            u_new[i_node] = val_old + params.dt * (diff[i_int] + reaction);
        }

        // Advance time
        t += params.dt;

        // 5) Boundary condition at new time
        u_new[0]             = params.A * std::sin(t) * std::sin(t);
        u_new[params.N-1]    = 0.0;

        // 6) Copy u_new -> u_old
        u_old = u_new;

        // 7) If we've reached the next output time, print a snapshot
        if (outputCount < params.P && t >= nextOutputTime - 1e-14) {
            print_snapshot(t, u_old, dx);
            outputCount++;
            nextOutputTime += outputInterval;
        }
    }

    return 0;
}
