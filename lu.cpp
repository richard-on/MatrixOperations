#include <cmath>
#include "lu.h"

LU::LU(Matrix a) {
    int exc = 0, row = 0;

    Matrix e = Matrix(a.length());
    for (int i = 0; i < a.length(); i++) {
        e(i, i) = 1;
    }


    for (int i = 0; i < a.length(); i++) {
        double uMax = 0;
        for (int r = i; r < a.length(); r++) {
            double Uii = a(r,i);
            int q = 0;
            while (q < i) {
                Uii -= a(r,q)*a(q,r);
                q++;
            }
            if(std::abs(Uii) > uMax) {
                uMax = std::abs(Uii);
                row = r;
            }
        }
        if (i != row) {
            exc++;
            for (int q = 0; q < a.length(); q++) {
                double temp = e(i, q);
                e(i, q) = e(row, q);
                e(row, q) = temp;
                temp = a(i,q);
                a(i,q) = a(row, q);
                a(row, q) = temp;
            }
        }

        int j = i;
        while (j < a.length()) {
            int q = 0;
            while (q < i) {
                a(i,j) -= a(i,q) * a(q,j);
                q++;
            }
            j++;
        }
        j = i + 1;
        while (j < a.length()) {
            int q = 0;
            while (q < i) {
                a(j,i) -= a(j,q) * a(q,i);
                q++;
            }
            a(j,i) = a(j,i) / a(i,i);
            j++;
        }
    }

    this->lu = a;
    this->p = e;

}

Vector LU::solve(const Vector& b) {
    Vector x(lu.length());
    Vector y(lu.length());

    Vector pb = p * b;
    for (int i = 0; i < lu.length(); i++) {
        y(i)=pb(i);
        int j = 0;
        while (j < i) {
            y(i) -= lu(i, j) * y(j);
            j++;
        }
    }
    for (int i = lu.length() - 1; i > -1; i--) {
        x(i)=y(i);
        for(int j = i+1; j < lu.length(); j++) {
            x(i) -= lu(i, j) * x(j);
        }
        x(i)= x(i) / lu(i, i);
    }

    return x;
}

const Matrix &LU::getLU() const {
    return lu;
}

const Matrix &LU::getP() const {
    return p;
}
