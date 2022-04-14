#ifndef MV_LDLT_H
#define MV_LDLT_H


#include "matrix.h"

class LDLT {
public:
    explicit LDLT(Matrix a);

    Vector solve(const Vector& b);

    const Matrix &getL() const;

    const Vector &getD() const;

private:
    Matrix l;
    Vector d;
};


#endif //MV_LDLT_H
