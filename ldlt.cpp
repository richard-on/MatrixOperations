#include "ldlt.h"

LDLT::LDLT(Matrix a) {
    Matrix lt = a;
    Vector d = Vector(3);

    /*for (int i = 0; i < a.length(); i++) {
        for (int j = 0; j < a.length(); j++) {
            if(i == r) {
                e(i, r) = 1;
            }
        }
    }*/
}

Vector LDLT::solve(const Vector &b) {
    return Vector();
}

const Matrix &LDLT::getL() const {
    return l;
}

const Vector &LDLT::getD() const {
    return d;
}
