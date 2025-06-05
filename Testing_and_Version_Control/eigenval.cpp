#include "eigenval.h"

// Internal routine to apply m to a, i.e., compute a <- m * a
void applyMatrix(const rmatrix<double>& m, rvector<double>& a) {
    int n = m.extent(1);
    rvector<double> q(n);
    q.fill(0.0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            q[i] += m[i][j]*a[j];
    a = q.copy();
}

// Internal routine to compute the expectation value of m in state a,
// i.e. <a|m|a>/<a|a>
double expectation(const rmatrix<double>& m, const rvector<double>& a) {
    int n = m.extent(1);
    double e = 0.0, z = 0.0;
    for (int i = 0; i < n; i++) {
        z += a[i]*a[i];
        for (int j = 0; j < n; j++) {
            e += a[i]*m[i][j]*a[j];
        }
    }
    if (z > 0)
        return e/z;
    else
        return e;
}

// Compute an approximation to the smallest eigenvalue e and
// corresponding eigenvector a of the matrix m.
// Not robust, may fail if large eigenvalues exist and will then throw.
void groundState(const rmatrix<double>& m,
                 double& e,
                 rvector<double>& a)
{
    int n = m.extent(1);
    // start with a uniform vector which will become the
    // (unnormalized) eigenvector
    rvector<double> eigenvector(n);
    eigenvector.fill(1.0);
    // find the largest diagonal element b, if positive, else 0
    double b = 0.0;
    for (int i = 0; i < n; i++)
        if (m[i][i] > b)
            b = m[i][i];
    // To find an approximation for the eigenvector with the lowest
    // energy, apply "m - b 1" many times to the eigenvector
    rmatrix<double> mshifted = m.copy();
    for (int i = 0; i < n; i++)
        mshifted[i][i] -= b;
    for (int p = 0; p < 20; p++)
        applyMatrix(mshifted, eigenvector);
    // the expectation value of this eigenvector is approximately the
    // ground state energy
    double eigenvalue = expectation(m, eigenvector);
    if (eigenvalue <= b) { // if this is not satisfied, the method failed
        e = eigenvalue;
        a = eigenvector;
    } else {
        std::cerr << eigenvalue << eigenvector;
        throw "Method failed";
    }
}
