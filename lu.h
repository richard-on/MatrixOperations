#ifndef MV_LU_H
#define MV_LU_H

#include "matrix.h"
#include "vector.h"

class lu {
public:
    explicit lu(Matrix a);

    Vector solve(const Vector& b);

    const Matrix &getData() const;

    const Matrix &getP() const;

private:
    Matrix e;
    Matrix data;
};


#endif //MV_LU_H
