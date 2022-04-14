#ifndef MV_LU_H
#define MV_LU_H

#include "matrix.h"
#include "vector.h"

class LU {
public:
    explicit LU(Matrix a);

    Vector solve(const Vector& b);

    const Matrix &getLU() const;

    const Matrix &getP() const;

private:
    Matrix lu;
    Matrix p;
};


#endif //MV_LU_H
