// output_c2.cpp
#include "output.h"
#include <fstream>

#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
using namespace Catch;

TEST_CASE("writeText test")
{
    // create file:
    rvector<double> a(3);
    a = 1,2,3;
    writeText("testoutputarr.txt", a);
    // read back:
    std::ifstream in("testoutputarr.txt");
    std::string s[3];
    in >> s[0] >> s[1] >> s[2];
    // check
    REQUIRE(s[0]=="1");
    REQUIRE(s[1]=="2");
    REQUIRE(s[2]=="3");
}

TEST_CASE("writeBinary test") {
    // Create rvector with values
    rvector<double> a(3);
    a = 1, 2, 3;

    // Write to binary file
    writeBinary("testoutputarr.bin", a);

    // Read back:
    rvector<double> b(3);
    std::ifstream in("testoutputarr.bin", std::ios::binary);
    in.read(reinterpret_cast<char*>(b.data()), b.size() * sizeof(double));

    // Check if the data was correctly written and read
    REQUIRE(b[0] == Approx(1.0));
    REQUIRE(b[1] == Approx(2.0));
    REQUIRE(b[2] == Approx(3.0));
}
