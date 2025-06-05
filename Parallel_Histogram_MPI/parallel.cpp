#include <mpi.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>

// Helper function for log base 'base'.
double logBase(double x, double base)
{
    return std::log(x) / std::log(base);
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // -----------------------------------------------------
    // 1) Parse command-line arguments (only on rank 0, then broadcast)
    // -----------------------------------------------------
    double base = 1.1;            // default log base
    std::string filename = "morestepnumbers_1.7GB.dat";
    int batchSize = 100000;       // default batch size (not strictly used here)
    if (rank == 0)
    {
        if (argc > 1) base     = std::stod(argv[1]);
        if (argc > 2) filename = argv[2];
        if (argc > 3) batchSize= std::stoi(argv[3]);
    }

    MPI_Bcast(&base,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&batchSize, 1, MPI_INT,    0, MPI_COMM_WORLD);
    // (We broadcast the batchSize as well, though we won't actually chunk the file here.)

    // -----------------------------------------------------
    // We'll do a 2-pass approach for min/max, then histogram
    // -----------------------------------------------------
    const int BINS = 30;
    double globalMin =  1e30;
    double globalMax = -1e30;

    // -----------------------------------------------------
    // PASS 1: Each rank reads its share of lines
    //         and finds localMin, localMax in log domain.
    // -----------------------------------------------------
    double localMin =  1e30;
    double localMax = -1e30;

    {
        std::ifstream fin(filename);
        if (!fin.good())
        {
            if (rank == 0)
            {
                std::cerr << "Error opening " << filename << std::endl;
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        long long lineIndex = 0; // Global line index
        double val;
        while (fin >> val)
        {
            // Only process lines that "belong" to this rank
            if ( (lineIndex % size) == rank )
            {
                double lx = logBase(val, base);
                if (lx < localMin) localMin = lx;
                if (lx > localMax) localMax = lx;
            }
            lineIndex++;
        }
        fin.close();
    }

    // Combine all localMin, localMax to find global min, max
    MPI_Allreduce(&localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&localMax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // -----------------------------------------------------
    // PASS 2: Build the local histogram, then reduce
    // -----------------------------------------------------
    std::vector<int> localHist(BINS, 0);
    long long localCount = 0;

    {
        std::ifstream fin(filename);
        if (!fin.good())
        {
            // This should not happen if pass 1 succeeded,
            // but let's be consistent.
            if (rank == 0)
            {
                std::cerr << "Error opening " << filename << std::endl;
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        long long lineIndex = 0;
        double val;
        while (fin >> val)
        {
            if ( (lineIndex % size) == rank )
            {
                double lx = logBase(val, base);
                int bin = int( (lx - globalMin) / (globalMax - globalMin) * BINS );
                if (bin < 0)       bin = 0;
                if (bin >= BINS)   bin = BINS - 1;
                localHist[bin]++;
                localCount++;
            }
            lineIndex++;
        }
        fin.close();
    }

    // Gather partial histograms and total count on rank 0
    std::vector<int> totalHist(BINS, 0);
    MPI_Reduce(localHist.data(), totalHist.data(), BINS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    long long globalCount = 0;
    MPI_Reduce(&localCount, &globalCount, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // -----------------------------------------------------
    // 4) Print the normalized histogram from rank 0
    // -----------------------------------------------------
    if (rank == 0)
    {
        for (int i = 0; i < BINS; i++)
        {
            // Bin i's left edge in the *log* domain
            double leftEdge = globalMin + (globalMax - globalMin) * (double(i) / BINS);
            double fraction = 0.0;
            if (globalCount > 0)
                fraction = (double)totalHist[i] / (double)globalCount;

            std::cout << leftEdge << " " << fraction << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
