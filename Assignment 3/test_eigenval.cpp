#define CATCH_CONFIG_MAIN
#include <catch2/catch_all.hpp>
#include "eigenval.h"
#include "rarray.hpp" // Required for rmatrix and rvector

using namespace Catch;

TEST_CASE("Testing groundState function", "[eigenval]") {
    SECTION("Ground state energy of {{-1,1},{2,0}} is -2") {
        rmatrix<double> mat(2, 2);
        mat[0][0] = -1; mat[0][1] = 1;
        mat[1][0] = 2;  mat[1][1] = 0;

        double e;
        rvector<double> a(2);
        groundState(mat, e, a);

        REQUIRE(e == Approx(-2.0));
    }

    SECTION("Ground state energy of identity matrix is 1") {
        rmatrix<double> mat(2, 2);
        mat.fill(0);
        mat[0][0] = 1;
        mat[1][1] = 1;

        double e;
        rvector<double> a(2);
        groundState(mat, e, a);

        REQUIRE(e == Approx(1.0));
    }

    SECTION("Ground state energy of diagonal matrix is its smallest element") {
        rmatrix<double> mat(3, 3);
        mat.fill(0);
        mat[0][0] = 3;
        mat[1][1] = 1;
        mat[2][2] = 2;

        double e;
        rvector<double> a(3);
        groundState(mat, e, a);

        REQUIRE(e == Approx(1.0));
    }

    SECTION("Function fails for matrix {{1,1},{2,0}}") {
        rmatrix<double> mat(2, 2);
        mat[0][0] = 1; mat[0][1] = 1;
        mat[1][0] = 2; mat[1][1] = 0;

        double e;
        rvector<double> a(2);

        REQUIRE_THROWS_AS(groundState(mat, e, a), const char*);
    }
}
