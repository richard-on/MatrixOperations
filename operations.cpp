#include <stdexcept>
#include <iostream>

#include "operations.h"

Vector operations::solveGauss(const Matrix& A, const Vector& b) {
    Vector x(b.length());

    auto** res = new double*[A.length()];
    for (int i = 0; i < A.length(); i++) {
        res[i] = new double[A.length()+ 1]{0};
        for(int j = 0; j < A.length() + 1; j++) {
            if(j == A.length()) {
                res[i][j] = b(i);
            }
            else {
                res[i][j] = A(i, j);
            }
        }
    }

    for (int i = 0; i < A.length() - 1; i++) {
        for (int j = i + 1; j < A.length(); j++) {
            double d = res[j][i]/res[i][i];
            for(int k = 0; k < A.length() + 1; k++) {
                res[j][k] = res[j][k]-d*res[i][k];
            }
        }
    }

    for(int i = A.length() - 1; i >= 0; i--) {
        x(i) = res[i][A.length()];
        for(int j = i + 1; j < A.length(); j++) {
            if(i != j) {
                x(i) = x(i) - (res[i][j] * x(j));
            }
        }
        x(i) = x(i) / res[i][i];
    }

    return x;
}

Vector operations::solveGaussCol(const Matrix& A, const Vector& b) {
    Vector x(b.length());

    double aa, bb;

    auto** res = new double*[A.length()];
    for (int i = 0; i < A.length(); i++) {
        res[i] = new double[A.length() + 1]{0};
        for(int j = 0; j < A.length() + 1; j++) {
            if(j == A.length()) {
                res[i][j] = b(i);
            }
            else {
                res[i][j] = A(i, j);
            }
        }
    }

    for (int i = 0, k = 0; k < A.length(); k++) {
        aa = abs(res[k][k]);
        i = k;
        for(int m = k + 1; m < A.length(); m++) {
            if (abs(res[m][k]) > aa) {
                i = m;
                aa = abs(res[m][k]);
            }
        }

        if (aa == 0) {
            throw std::invalid_argument("system doesn't have any solutions");
        }

        if (i != k) {
            for (int j = k; j < A.length() + 1; j++) {
                bb = res[k][j];
                res[k][j] = res[i][j];
                res[i][j] = bb;
            }
        }

        aa = res[k][k];
        res[k][k] = 1;
        for (int j = k + 1; j < A.length() + 1; j++) {
            res[k][j] = res[k][j] / aa;
        }
        for (i = k + 1; i < A.length(); i++) {
            bb = res[i][k];
            res[i][k] = 0;
            if (bb != 0) {
                for (int j = k + 1; j < A.length() + 1; j++) {
                    res[i][j] = res[i][j] - bb * res[k][j];
                }
            }
        }
    }

    for(int i = A.length() - 1; i >= 0; i--) {
        x(i) = res[i][A.length()];
        for(int j = i + 1; j < A.length(); j++) {
            if(i != j) {
                x(i) = x(i) - (res[i][j] * x(j));
            }
        }
        x(i) = x(i) / res[i][i];
    }

    return x;
}

Vector operations::solveSOR(const Matrix &A, const Vector &b, double param, double eps) {
    Vector x(b.length());

    int step = 0;
    Vector m = A * x;
    double res = (m - b).norm();

    while(res > eps) {
        for(int i = 0; i < A.length(); i++) {
            double norm = 0;
            for(int j = 0; j < A.length(); j++) {
                if(j != i) {
                    norm += A(i,j)*x(j);
                }
            }
            x(i) = (1 - param) * x(i) + (param / A(i,i)) * (b(i) - norm);
        }
        m = A * x;
        res = (m - b).norm();
        step++;
        //std::cout << "Step: " << step << " Res: " << res << std::endl;
    }

    return x;
}
