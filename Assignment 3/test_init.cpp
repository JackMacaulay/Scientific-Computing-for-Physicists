#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "init.h"
using namespace Catch;

TEST_CASE("Testing initMatrix function", "[init]") {
    SECTION("n is zero") {
        rmatrix<double> result = initMatrix(0);
        REQUIRE(result.size() == 0);  
    }

    SECTION("n is one") {
        rmatrix<double> result = initMatrix(1);
        REQUIRE(result.size() == 1);
        REQUIRE(result[0][0] == Approx(-97.0 / 2500.0));  
    }

    SECTION("n is the cube of an odd number") {
        int n = 27; 
        rmatrix<double> result = initMatrix(n);
        REQUIRE(result.size() == n*n);
        
    }

    SECTION("n is not the cube of an odd number") {
        int n = 8; 
        REQUIRE_THROWS_AS(initMatrix(n), std::invalid_argument);
    }
}
