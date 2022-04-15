#ifndef MV_LDLT_H
#define MV_LDLT_H


#include "matrix.h"

class LDLT {
public:
    explicit LDLT(const Matrix& a);

    Vector solve(const Vector& b);

    const Matrix &getLLT() const;

    const Vector &getD() const;

private:
    Matrix LLT;
    Vector D;
};


#endif //MV_LDLT_H
