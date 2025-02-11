#define CATCH_CONFIG_MAIN
#include <Catch2/catch_all.hpp>
#include "init.h"
#include "eigenval.h"
#include "output.h"

TEST_CASE("Integrated test of the entire code", "[integration]") {
    int n = 27; // 3^3
    auto matrix = init(n);
    REQUIRE(matrix.size() == n);

    double ground_state_energy = eigenval(matrix);
    // Add checks for ground state energy

    std::vector<double> data = {ground_state_energy};
    std::string filename = "integration_test_output.bin";
    writeBinary(data, filename);

    // Verify the file contents
    std::ifstream infile(filename, std::ios::binary);
    REQUIRE(infile.is_open());

    std::vector<double> read_data(data.size());
    infile.read(reinterpret_cast<char*>(read_data.data()), read_data.size() * sizeof(double));
    REQUIRE(infile.gcount() == static_cast<std::streamsize>(data.size() * sizeof(double)));

    REQUIRE(read_data == data);

    // Clean up
    infile.close();
    std::remove(filename.c_str());
}
