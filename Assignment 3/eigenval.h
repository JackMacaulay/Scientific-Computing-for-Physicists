#ifndef EIGENVALH
#define EIGENVALH

#include <rarray>

// Compute an approximation to the smallest eigenvalue e and
// corresponding unnormalized eigenvector a of the matrix m.
// Not robust, may fail if large eigenvalues exist and will then throw.
void groundState(const rmatrix<double>& m, double& e, rvector<double>& a);

#endif
