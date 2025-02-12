#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "init.h"
#include "eigenval.h"
#include "output.h"
#include <fstream>
#include "rarray.hpp"  // Needed for rmatrix and rvector

using namespace Catch;

TEST_CASE("Integrated test of the entire code", "[integration]") {
    // Verify matrix is correct size
    int n = 27;
    rmatrix<double> matrix = initMatrix(n); 
    REQUIRE(matrix.size() == n*n);

    // Compute ground state energy
    double ground_state_energy;
    rvector<double> eigenvector(n);
    groundState(matrix, ground_state_energy, eigenvector); 

    
    REQUIRE(ground_state_energy < 0.0); // Ground state energy should be negative

    // Prepare data for writing
    rvector<double> data(1);
    data[0] = ground_state_energy;
    std::string filename = "integration_test_output.bin";

    // Write binary data
    writeBinary(filename, data);

    // Verify file contents
    std::ifstream infile(filename, std::ios::binary);
    REQUIRE(infile.is_open());

    rvector<double> read_data(1);
    infile.read(reinterpret_cast<char*>(read_data.data()), read_data.size() * sizeof(double));

    REQUIRE(infile.gcount() == static_cast<std::streamsize>(read_data.size() * sizeof(double)));
    REQUIRE(read_data[0] == Approx(ground_state_energy));

    // Clean up
    infile.close();
    std::remove(filename.c_str());
}
