#include <cmath>
#include "lu.h"
#include "operations.h"

lu::lu(Matrix a) {
    int exc = 0, row = 0;

    Matrix e = Matrix(a.length());
    for (int i = 0; i < a.length(); i++) {
        for (int r = 0; r < a.length(); r++) {
            if(i == r) {
                e(i, r) = 1;
            }
        }
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

    this->data = a;
    this->e = e;

}

Vector lu::solve(const Vector& b) {
    Vector x(data.length());
    Vector y(data.length());

    Vector pb = e * b;
    for (int i = 0; i < data.length(); i++) {
        y(i)=pb(i);
        int j = 0;
        while (j < i) {
            y(i) -= data(i, j)*y(j);
            j++;
        }
    }
    for (int i = data.length() - 1; i > -1; i--) {
        x(i)=y(i);
        int j = i+1;
        while (j < data.length()) {
            x(i) -= data(i, j)*x(j);
            j++;
        }
        x(i)=x(i)/data(i, i);
    }

    return x;
}

const Matrix &lu::getData() const {
    return data;
}

const Matrix &lu::getP() const {
    return e;
}
