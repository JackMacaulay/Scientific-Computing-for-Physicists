// hydrogen.cpp
#include <iostream>
#include <rarray>
#include "eigenval.h"
#include "output.h"
#include "init.h"
int main() {
    const int n = 4913;
    rmatrix<double> m = initMatrix(n);
    rvector<double> a;
    double e;    
    groundState(m, e, a);
    std::cout<<"Ground state energy="<<e<<"\n";
    writeText("data.txt", a);
    writeBinary("data.bin", a);
}
