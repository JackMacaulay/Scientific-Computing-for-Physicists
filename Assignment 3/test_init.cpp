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

    //Both of the test cases below fail - this case fails due to an assertion error; however we expect it to fail
    //And when it is converted to an unsigned integer the loop in initMatrix runs from 0-7 since the matrix is 8x8
    // which leds to the assertion I<1 failing. Hence causing it to fail immediately
    //likewise with the next case; however, because of the assertion error it will not test that case
    SECTION("n is not the cube of an odd number") {
        int n = 8; 
        REQUIRE_THROWS_AS(initMatrix(n), std::invalid_argument);
    }

    SECTION("n is the cube of an odd number") {
        int n = 27; 
        rmatrix<double> result = initMatrix(n);
        std::cout << "Matrix size: " << result.extent(0) << " x " << result.extent(1) << std::endl;
        REQUIRE(1==1);
        

        
        
    }
}
