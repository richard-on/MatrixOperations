#include <iostream>
#include <cmath>
#include "ldlt.h"

LDLT::LDLT(const Matrix& m) {

    Matrix a = Matrix(m);
    Vector d = Vector(a.length());

    for (int i = 0; i < a.length() - 1; i++) {
        for (int j = i + 1; j < a.length(); j++) {
            double div = a(j, i) / a(i, i);
            for(int k = 0; k < a.length(); k++) {
                a(j, k) -= div * a(i, k);
            }
        }
    }

    int sign;
    for(int i = 0; i < a.length(); i++) {
        double div = std::sqrt(std::abs(a(i, i)));
        if(a(i, i) >= 0) {
            sign = 1;
        } else{
            sign = -1;
        }
        d(i) = sign;
        a(i,i) /= div;
        for(int j = i + 1; j < a.length(); j++) {
            a(i,j) /= div;
            a(j,i) = a(i,j);
            a(i,j) *= sign;
        }
    }

    this->LLT = a;
    this->D = d;
}

Vector LDLT::solve(const Vector &b) {
    Vector x(LLT.length());
    Vector y(LLT.length());

    double sum;
    for(int i = 0; i < LLT.length(); i++) {
        sum = 0;
        for(int j = 0; j <= i; j++) {
            sum += LLT(i, j) * y(j);
        }
        y(i) = (b(i) - sum) / LLT(i,i);
    }

    for(int i = LLT.length() - 1; i >= 0; i--) {
        sum = 0;
        for(int j = LLT.length() - 1; j >= i; j--) {
            sum += LLT(i, j) * x(j);
        }
        x(i) = (y(i) - sum) / LLT(i,i);
    }

    return x;
}

const Matrix &LDLT::getLLT() const {
    return LLT;
}

const Vector &LDLT::getD() const {
    return D;
}
