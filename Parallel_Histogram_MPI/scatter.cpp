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
    // 1) Parse command-line arguments (only on root).
    // -----------------------------------------------------
    double base = 1.1;         // default log base
    std::string filename = "morestepnumbers.dat";
    int batchSize = 100000;    // default batch size
    if (rank == 0)
    {
        if (argc > 1) base     = std::stod(argv[1]);
        if (argc > 2) filename = argv[2];
        if (argc > 3) batchSize= std::stoi(argv[3]);
    }

    // Broadcast the numeric arguments to all ranks.
    MPI_Bcast(&base,     1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&batchSize,1, MPI_INT,    0, MPI_COMM_WORLD);
    // (Filename is only opened on rank 0 below.)

    // We'll store a local histogram on each rank and later reduce them.
    const int BINS = 30;
    std::vector<int> localHist(BINS, 0);
    long long localCount = 0;  // number of data points processed on this rank

    // -----------------------------------------------------
    // First pass: find global min & max of log-values
    // -----------------------------------------------------
    double globalMin =  1e30;
    double globalMax = -1e30;

    {
        std::ifstream fin;
        if (rank == 0)
        {
            fin.open(filename);
            if (!fin)
            {
                std::cerr << "Error opening " << filename << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        while (true)
        {
            // Root reads up to batchSize numbers from file
            std::vector<double> batch;
            int count = 0;
            if (rank == 0)
            {
                batch.reserve(batchSize);
                double val;
                while (count < batchSize && (fin >> val))
                {
                    batch.push_back(val);
                    count++;
                }
            }

            // Number actually read in this batch:
            MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

            // If no more data, break
            if (count == 0) {
                break;
            }

            // Because 'count' may not be divisible evenly by 'size',
            // we use Scatterv for uneven distribution.
            std::vector<int> sendCounts(size), displs(size);
            {
                int q  = count / size;
                int r  = count % size;
                for (int i = 0; i < size; i++)
                    sendCounts[i] = q + (i < r ? 1 : 0);

                displs[0] = 0;
                for (int i = 1; i < size; i++)
                    displs[i] = displs[i - 1] + sendCounts[i - 1];
            }

            int localN = sendCounts[rank];
            std::vector<double> localData(localN, 0.0);

            MPI_Scatterv(batch.data(), sendCounts.data(), displs.data(),
                         MPI_DOUBLE,
                         localData.data(), localN, MPI_DOUBLE,
                         0, MPI_COMM_WORLD);

            // Compute local min/max in the log domain
            double myMin =  1e30;
            double myMax = -1e30;
            for (double x : localData)
            {
                double lx = logBase(x, base);
                if (lx < myMin) myMin = lx;
                if (lx > myMax) myMax = lx;
            }

            // Combine to get global min & max
            double tmp;
            MPI_Allreduce(&myMin, &tmp, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
            if (tmp < globalMin) globalMin = tmp;

            MPI_Allreduce(&myMax, &tmp, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            if (tmp > globalMax) globalMax = tmp;
        }

        if (rank == 0) fin.close();
    }

    // -----------------------------------------------------
    // Second pass: actually build the histograms
    // -----------------------------------------------------
    {
        std::ifstream fin;
        if (rank == 0)
        {
            fin.open(filename);
            if (!fin)
            {
                std::cerr << "Error opening " << filename << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        while (true)
        {
            std::vector<double> batch;
            int count = 0;
            if (rank == 0)
            {
                batch.reserve(batchSize);
                double val;
                while (count < batchSize && (fin >> val))
                {
                    batch.push_back(val);
                    count++;
                }
            }

            MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);

            if (count == 0) {
                break;
            }

            std::vector<int> sendCounts(size), displs(size);
            {
                int q  = count / size;
                int r  = count % size;
                for (int i = 0; i < size; i++)
                    sendCounts[i] = q + (i < r ? 1 : 0);

                displs[0] = 0;
                for (int i = 1; i < size; i++)
                    displs[i] = displs[i - 1] + sendCounts[i - 1];
            }

            int localN = sendCounts[rank];
            std::vector<double> localData(localN, 0.0);

            MPI_Scatterv(batch.data(), sendCounts.data(), displs.data(),
                         MPI_DOUBLE,
                         localData.data(), localN, MPI_DOUBLE,
                         0, MPI_COMM_WORLD);

            // Bin each local point
            for (double x : localData)
            {
                double lx = logBase(x, base);
                int bin = int((lx - globalMin) / (globalMax - globalMin) * BINS);
                if (bin < 0)      bin = 0;
                if (bin >= BINS)  bin = BINS - 1;
                localHist[bin]++;
                localCount++;
            }
        }

        if (rank == 0) fin.close();
    }

    // -----------------------------------------------------
    // 4) Reduce all partial histograms to rank 0 & print
    // -----------------------------------------------------
    std::vector<int> totalHist(BINS, 0);
    MPI_Reduce(localHist.data(), totalHist.data(),
               BINS, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    long long globalCount = 0;
    MPI_Reduce(&localCount, &globalCount,
               1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Print the normalized histogram from rank 0
    if (rank == 0)
    {
        for (int i = 0; i < BINS; i++)
        {
            // Left edge of this bin in the log domain:
            double leftEdge = globalMin + (globalMax - globalMin) * (double(i) / BINS);
            double fraction = double(totalHist[i]) / double(globalCount);
            std::cout << leftEdge << " " << fraction << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
