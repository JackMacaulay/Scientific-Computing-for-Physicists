#define CATCH_CONFIG_MAIN
#include <Catch2/catch_all.hpp>
#include "init.h"

TEST_CASE("Testing init function", "[init]") {
    SECTION("n is zero") {
        auto result = init(0);
        REQUIRE(result.empty());
    }

    SECTION("n is one") {
        auto result = init(1);
        REQUIRE(result.size() == 1);
        REQUIRE(result(0, 0) == Approx(-97.0 / 2500.0));
    }

    SECTION("n is the cube of an odd number") {
        int n = 27; // 3^3
        auto result = init(n);
        REQUIRE(result.size() == n);
        // Additional checks can be added here
    }

    SECTION("n is not the cube of an odd number") {
        int n = 8; // 2^3
        REQUIRE_THROWS_AS(init(n), std::invalid_argument);
    }
}
