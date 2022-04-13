#ifndef MV_OPERATIONS_H
#define MV_OPERATIONS_H


#include "vector.h"
#include "matrix.h"

class operations {
public:
    static Vector solveGauss(const Matrix& A, const Vector& b);
    static Vector solveGaussCol(const Matrix& A, const Vector& b);
    static Vector solveSOR(const Matrix& A, const Vector& b, double param, double eps = 1e-8);
};


#endif //MV_OPERATIONS_H
