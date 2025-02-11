// init.cpp
#include "init.h"
#include <cassert>
#include <cmath>

double hammat(unsigned I, unsigned J, double L, unsigned r)
{
    const unsigned n = 2*r+1;
    const unsigned nnn = n*n*n;    
    assert(I < nnn);
    assert(J < nnn);
    const unsigned i1 = I%n;
    const unsigned j1 = (I/n)%n;
    const unsigned k1 = (I/n)/n;
    const unsigned i2 = J%n;
    const unsigned j2 = (J/n)%n;
    const unsigned k2 = (J/n)/n;
    const double dx = L/n;
    const double xc1 = (i1-r)*dx;
    const double yc1 = (j1-r)*dx;
    const double zc1 = (k1-r)*dx;
    auto V = [](double x, double y, double z) { return -1.0/sqrt(x*x+y*y+z*z); };
    if (i1==i2 && j1==j2 && k1==k2) {
        const double Vijk = (V(xc1-dx/2,yc1,zc1)+V(xc1+dx/2,yc1,zc1)
                             +V(xc1,yc1-dx/2,zc1)+V(xc1,yc1+dx/2,zc1)
                             +V(xc1,yc1,zc1-dx/2)+V(xc1,yc1,zc1+dx/2))/6;
        const double t = -1/(2*dx*dx);
        return Vijk - 6*t;
    } else if ( (i1==i2 && j1==j2 && std::abs(static_cast<int>(k1 - k2)) == 1)
                || (i1==i2 && std::abs(static_cast<int>(j1 - j2)) == 1 && k1==k2)
                || (std::abs(static_cast<int>(i1 - i2)) == 1 && j1==j2 && k1==k2) ) {
        const double t = -1/(2*dx*dx);
        return t;
    } else {
        return 0;
    }
}

double hamiltonian(unsigned I, unsigned J, unsigned n, double L)
{
    unsigned r = (cbrt(n)-1)/2;
    while ((2*r+1)*(2*r+1)*(2*r+1) > n) r--;
    return hammat(I,J,L,r);
}

double H(int i,int j,int n) {
    const double L=50.0;
    return hamiltonian(i,j,n,L);
}


rmatrix<double> initMatrix(int n)
{
    rmatrix<double> m(n,n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            m[i][j] = H(i,j,n);
        }
    }
    return m;
}
