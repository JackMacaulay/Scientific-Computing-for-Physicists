
#include <iostream>
#include <iomanip>  
#include <cmath>
#include <cstdlib>
#include <complex>
#include <rarray>   
#include <fftw3.h>
#include <string>

// Structure to store simulation parameters.
struct SimulationParams {
    int P;       // Number of snapshots
    double L;    // Domain length
    double A;    // Amplitude in the initial condition
    int N;       // Number of grid points
    double T;    // Final time
    double dt;   // Time step
};

// Parse command-line arguments: expected: P L A N T dt
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

// Reaction substep: exact update for logistic growth.
inline double reaction_update(double u, double dt) {
    if(u <= 0.0) return 0.0;
    if(u >= 1.0) return 1.0;
    double exp_dt = std::exp(dt);
    return (u * exp_dt) / ((1.0 - u) + u * exp_dt);
}

// Print a snapshot: each line prints time, x, and u(x)
void print_snapshot(double t, const rarray<double,1>& u, double dx) {
    int N = u.size();
    for (int i = 0; i < N; i++) {
        std::cout << std::fixed << std::setprecision(6)
                  << t << " " << (i * dx) << " " << u[i] << "\n";
    }
}

int main(int argc, char* argv[]) {
    // Get simulation parameters.
    SimulationParams params = parse_arguments(argc, argv);
    int N = params.N;
    double L = params.L;
    double dt = params.dt;
    double T = params.T;
    int P = params.P;
    double A = params.A;
    double dx = L / N;
    
    // Create solution array u (size N) as an rvector.
    rarray<double,1> u(N);
    // Initial condition: u(x,0) = A * [ sin(π*x/L) ]^100.
    for (int j = 0; j < N; ++j) {
        double x = j * dx;
        u[j] = A * std::pow(std::sin(M_PI * x / L), 100);
    }
    
    // Set up FFTW for the diffusion substep.
    // Allocate FFT input and output arrays.
    double* fft_in = fftw_alloc_real(N);
    fftw_complex* fft_out = fftw_alloc_complex(N/2 + 1);
    fftw_plan plan_forward = fftw_plan_dft_r2c_1d(N, fft_in,
                                    reinterpret_cast<fftw_complex*>(fft_out),
                                    FFTW_ESTIMATE);
    fftw_plan plan_inverse = fftw_plan_dft_c2r_1d(N,
                                    reinterpret_cast<fftw_complex*>(fft_out),
                                    fft_in, FFTW_ESTIMATE);
    
    // Set up snapshot output.
    double t_current = 0.0;
    int snapshotCount = 0;
    double outputInterval = (P > 1) ? T / (P - 1.0) : T;
    double nextOutputTime = 0.0;
    
    // Print initial snapshot at t=0.
    if (P > 0) {
        print_snapshot(t_current, u, dx);
        snapshotCount++;
        nextOutputTime += outputInterval;
    }
    
    // Main time-stepping loop.
    const double EPS = 1e-14;
    while (t_current < T - EPS) {
        // Adjust dt for the last step if necessary.
        double dt_local = dt;
        if (t_current + dt_local > T)
            dt_local = T - t_current;
        
        // --- Reaction Substep ---
        // Update u at every grid point (including boundaries if needed).
        // For periodic boundaries, u(0) and u(N-1) must eventually satisfy u(0)=u(N-1).
        // Here we update all points using the logistic update.
        for (int j = 0; j < N; ++j) {
            u[j] = reaction_update(u[j], dt_local);
        }
        
        // --- Diffusion Substep ---
        // Copy u into fft_in.
        for (int j = 0; j < N; ++j)
            fft_in[j] = u[j];
        
        // Forward FFT.
        fftw_execute(plan_forward);
        // Update each Fourier mode: u_hat(k,t+dt) = u_hat(k,t) * exp(-k^2 * dt_local)
        // Here, k = 2π*k_index/L for k_index = 0,1,...,N/2.
        for (int k = 0; k <= N/2; ++k) {
            double k_val = 2.0 * M_PI * k / L;
            double factor = std::exp(- (k_val * k_val) * dt_local);
            fft_out[k][0] *= factor;  // real part
            fft_out[k][1] *= factor;  // imaginary part
        }
        // Inverse FFT.
        fftw_execute(plan_inverse);
        // Normalize the inverse FFT (FFTW does not include normalization).
        for (int j = 0; j < N; ++j) {
            u[j] = fft_in[j] / N;
        }
        
        // For periodic BCs, ensure that the boundaries match.
        // Here we simply set u[0] = u[N-1] as the average of the two.
        double avgBoundary = 0.5 * (u[0] + u[N-1]);
        u[0] = avgBoundary;
        u[N-1] = avgBoundary;
        
        // Advance simulation time.
        t_current += dt_local;
        
        // Print snapshot if we've reached (or exceeded) the next output time.
        if (snapshotCount < P && t_current >= nextOutputTime - EPS) {
            print_snapshot(t_current, u, dx);
            snapshotCount++;
            nextOutputTime += outputInterval;
        }
    }
    
    // Ensure final snapshot at t=T is printed.
    if (snapshotCount < P) {
        t_current = T;
        print_snapshot(t_current, u, dx);
    }
    
    // Clean up FFTW resources.
    fftw_destroy_plan(plan_forward);
    fftw_destroy_plan(plan_inverse);
    fftw_free(fft_in);
    fftw_free(fft_out);
    
    return 0;
}
