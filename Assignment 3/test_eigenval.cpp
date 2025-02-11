#define CATCH_CONFIG_MAIN
#include <Catch2/catch_all.hpp>
#include "eigenval.h"

TEST_CASE("Testing eigenval function", "[eigenval]") {
    SECTION("Ground state energy of {{-1,1},{2,0}} is -2") {
        Eigen::Matrix2d mat;
        mat << -1, 1,
               2, 0;
        REQUIRE(eigenval(mat) == Approx(-2.0));
    }

    SECTION("Ground state energy of identity matrix is 1") {
        Eigen::Matrix2d mat = Eigen::Matrix2d::Identity();
        REQUIRE(eigenval(mat) == Approx(1.0));
    }

    SECTION("Ground state energy of diagonal matrix is its smallest element") {
        Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
        mat.diagonal() << 3, 1, 2;
        REQUIRE(eigenval(mat) == Approx(1.0));
    }

    SECTION("Function fails for matrix {{1,1},{2,0}}") {
        Eigen::Matrix2d mat;
        mat << 1, 1,
               2, 0;
        REQUIRE_THROWS_AS(eigenval(mat), std::runtime_error);
    }
}
