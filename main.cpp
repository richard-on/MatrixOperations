#include <iostream>
#include <random>
#include <iomanip>
#include <cmath>

#define VAR 9
#define LENGTH 3

double** multiplyMatrix(auto& a, auto& b) {
    auto** res = new double*[LENGTH];
    for(int i = 0; i < LENGTH; i++) {
        res[i] = new double[LENGTH]{0};
    }

    for(int i = 0; i < LENGTH; i++) {
        for(int j = 0; j < LENGTH; j++) {
            for(int k = 0; k < LENGTH; k++) {
                res[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return res;
}

void printMatrix(auto &m, int r, int c) {
    for (int i = 0; i < r; i++) {
        if(i == 0) {
            std::cout << "/";
        }
        else if(i == r - 1) {
            std::cout << "\\";
        }
        else{
            std::cout << "|";
        }

        for (int j = 0; j < c; j++) {

            std::cout << std::setw(8) << std::setprecision(3) << m[i][j] << " ";
        }
        if(i == 0) {
            std::cout << "\\";
        }
        else if(i == r - 1) {
            std::cout << "/";
        }
        else{
            std::cout << "|";
        }
        std::cout << std::endl;

    }
    std::cout << std::endl;
}

void printVector(double** vec) {
    for (int i = 0; i < LENGTH; i++) {
        if(i == 0) {
            std::cout << "/";
        }
        else if(i == LENGTH - 1) {
            std::cout << "\\";
        }
        else{
            std::cout << "|";
        }
        std::cout << std::setw(6) << std::setprecision(3) << vec[i][0];
        if(i == 0) {
            std::cout << "\\";
        }
        else if(i == LENGTH - 1) {
            std::cout << "/";
        }
        else{
            std::cout << "|";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

double** inverseMatrix(double** a) {
    auto** res = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        res[i] = new double[LENGTH]{0};
    }


    auto** m = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        m[i] = new double[LENGTH * 2]{0};
        for (int j = 0; j < 2 * LENGTH; j++) {
            if(j < LENGTH) {
                m[i][j] = a[i][j];
            }
            if (j == (i + LENGTH))
                m[i][j] = 1;
        }
    }

    for (int i = LENGTH-1; i > 0; i--) {
        if (m[i - 1][0] < m[i][0]) {
            double* temp = m[i];
            m[i] = m[i - 1];
            m[i - 1] = temp;
        }
    }

    for (int i = 0; i < LENGTH; i++) {
        for (int j = 0; j < LENGTH; j++) {
            if (j != i) {
                double d = m[j][i] / m[i][i];
                for (int k = 0; k < 2 * LENGTH; k++) {
                    m[j][k] -= m[i][k] * d;
                }
            }
        }
    }

    for (int i = 0; i < LENGTH; i++) {
        double d = m[i][i];
        for (int j = 0; j < 2 * LENGTH; j++) {
            m[i][j] = m[i][j] / d;
        }
    }

    for (int i = 0; i < LENGTH; i++) {
        int k = 0;
        for (int j = LENGTH; j < 2 * LENGTH; j++) {
            res[i][k] = m[i][j];
            k++;
        }
    }

    delete[] m;
    return res;
}

double cubicNorm(double** m) {
    double norm = INT_MIN;
    double columnSum = 0;

    for (int i = 0; i < LENGTH; i++) {
        for (int j = 0; j < LENGTH; j++) {
            columnSum += std::abs(m[i][j]);
        }
        if(columnSum > norm) {
            norm = columnSum;
        }
        columnSum = 0;
    }

    return norm;
}

double** solveGaussCol(double** A, double** b) {
    double aa, bb;

    auto** x = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        x[i] = new double[1]{0};
    }

    auto** res = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        res[i] = new double[LENGTH + 1]{0};
        for(int j = 0; j < LENGTH + 1; j++) {
            if(j == LENGTH) {
                res[i][j] = b[i][0];
            }
            else {
                res[i][j] = A[i][j];
            }
        }
    }

    for (int i = 0, k = 0; k < LENGTH; k++) {
        aa = abs(res[k][k]);
        i = k;
        for(int m = k + 1; m < LENGTH; m++) {
            if (abs(res[m][k]) > aa) {
                i = m;
                aa = abs(res[m][k]);
            }
        }

        if (aa == 0) {
            std::cout<<"Система не имеет решений"<<std::endl;
        }

        if (i != k) {
            for (int j = k; j < LENGTH + 1; j++) {
                bb = res[k][j];
                res[k][j] = res[i][j];
                res[i][j] = bb;
            }
        }

        aa = res[k][k];
        res[k][k] = 1;
        for (int j = k + 1; j < LENGTH + 1; j++) {
            res[k][j] = res[k][j] / aa;
        }
        for (i = k + 1; i < LENGTH; i++) {
            bb = res[i][k];
            res[i][k] = 0;
            if (bb != 0) {
                for (int j = k + 1; j < LENGTH + 1; j++) {
                    res[i][j] = res[i][j] - bb * res[k][j];
                }
            }
        }
    }

    for(int i = LENGTH - 1; i >= 0; i--) {
        x[i][0] = res[i][LENGTH];
        for(int j = i + 1; j < LENGTH; j++) {
            if(i != j) {
                x[i][0] = x[i][0] - (res[i][j] * x[j][0]);
            }
        }
        x[i][0]= x[i][0] / res[i][i];
    }

    return x;
}

double** solveGauss(double** A, double** b) {
    auto** x = new double*[LENGTH];
    for(int i = 0; i < LENGTH; i++) {
        x[i] = new double[1]{0};
    }

    auto** res = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        res[i] = new double[LENGTH + 1]{0};
        for(int j = 0; j < LENGTH + 1; j++) {
            if(j == LENGTH) {
                res[i][j] = b[i][0];
            }
            else {
                res[i][j] = A[i][j];
            }
        }
    }

    for (int i = 0; i < LENGTH - 1; i++) {
        for (int j = i + 1; j < LENGTH; j++) {
            double d = res[j][i]/res[i][i];
            for(int k = 0; k < LENGTH + 1; k++) {
                res[j][k] = res[j][k]-d*res[i][k];
            }
        }
    }

    for(int i = LENGTH - 1; i >= 0; i--) {
        x[i][0] = res[i][LENGTH];
        for(int j = i + 1; j < LENGTH; j++) {
            if(i != j) {
                x[i][0] = x[i][0] - (res[i][j] * x[j][0]);
            }
        }
        x[i][0]= x[i][0] / res[i][i];
    }

    return x;
}

struct LU{
    double** lu;
    double** p;
};

double** solveLU(LU lu, double** b) {
    auto** x = new double*[LENGTH];
    for(int i = 0; i < LENGTH; i++) {
        x[i] = new double[1]{0};
    }
    auto** y = new double*[LENGTH];
    for(int i = 0; i < LENGTH; i++) {
        y[i] = new double[1]{0};
    }
    auto** pb = multiplyMatrix(lu.p, b);
    for (int q = 0; q < 1; q++) {
        for (int i = 0; i < LENGTH; i++) {
            y[i][q]=pb[i][q];
            int j = 0;
            while (j < i) {
                y[i][q] -= lu.lu[i][j]*y[j][q];
                j++;
            }
        }
        for (int i = LENGTH - 1; i > -1; i--) {
            x[i][q]=y[i][q];
            int j = i+1;
            while (j < LENGTH) {
                x[i][q] -= lu.lu[i][j]*x[j][q];
                j++;
            }
            x[i][q]=x[i][q]/lu.lu[i][i];
        }
    }
    delete[] pb;

    return x;
}

LU luFactorize(double** A) {
    LU lu{};
    int exc = 0, row = 0;
    auto** P = new double*[LENGTH];
    for (int i = 0; i < LENGTH; i++) {
        P[i] = new double[LENGTH]{0};
        for (int r = 0; r < LENGTH; r++) {
            if(i == r) {
                P[i][r] = 1;
            }
        }
    }

    for (int i = 0; i < LENGTH; i++) {
        double uMax = 0;
        for (int r = i; r < LENGTH; r++) {
            double Uii = A[r][i];
            int q = 0;
            while (q < i) {
                Uii -= A[r][q]*A[q][r];
                q++;
            }
            if(abs(Uii) > uMax) {
                uMax = abs(Uii);
                row = r;
            }
        }
        if (i != row) {
            exc++;
            for (int q = 0; q < LENGTH; q++) {
                double temp = P[i][q];
                P[i][q] = P[row][q];
                P[row][q] = temp;
                temp = A[i][q];
                A[i][q] = A[row][q];
                A[row][q] = temp;
            }
        }

        int j = i;
        while (j < LENGTH) { //Determine U across row i
            int q = 0;
            while (q < i) {
                A[i][j] -= A[i][q] * A[q][j];
                q++;
            }
            j++;
        }
        j = i+1;
        while (j < LENGTH) { //Determine L down column i
            int q = 0;
            while (q < i) {
                A[j][i] -= A[j][q] * A[q][i];
                q++;
            }
            A[j][i] = A[j][i] / A[i][i];
            j++;
        }
    }
    lu.lu = A;
    lu.p = P;

    return lu;
}

int main() {
    //int low_interval = static_cast<int>(pow(-2, VAR/4));
    //int high_interval = static_cast<int>(pow(2, VAR/4));


    /*std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_int_distribution<> distribution(-1, 1);*/

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution(-3, 3);

    auto** matrix = new double*[LENGTH];
    double rowSum = 0;
    for(int i = 0; i < LENGTH; i++) {
        matrix[i] = new double[LENGTH]{0};
        for(int j = LENGTH - 1; j >= 0; j--) {
            if(j > i) {
                matrix[i][j] = distribution(gen);
            }
            else if(j < i) {
                matrix[i][j] = matrix[j][i];
            }
            rowSum += abs(matrix[i][j]);

            if(j == 0){
                matrix[i][i] = rowSum + 1;
                rowSum = 0;
            }
        }
    }

    auto* y = new double*[LENGTH];
    for(int i = 0; i < LENGTH; i++) {
        y[i] = new double[1]{0};
        y[i][0] = distribution(gen);
    }

    std::cout << "Initial matrix: " << std::endl;
    printMatrix(matrix, LENGTH, LENGTH);

    std::cout << "Vector y: " << std::endl;
    printVector(y);

    double** b = multiplyMatrix(matrix, y);
    std::cout << "Vector b: " << std::endl;
    printVector(b);

    double** inverse = inverseMatrix(matrix);
    std::cout << "Inverse matrix: " << std::endl;
    printMatrix(inverse, LENGTH, LENGTH);

    double mc = cubicNorm(matrix);
    double mi = cubicNorm(inverse);
    double cond = mc * mi;
    std::cout << "Cond: " << cond << std::endl;

    double** gauss = solveGauss(matrix, b);
    std::cout << "Gauss solution: " << std::endl;
    printVector(gauss);

    double** gaussCol = solveGaussCol(matrix, b);
    std::cout << "Gauss Column solution: " << std::endl;
    printVector(gaussCol);

    LU lu = luFactorize(matrix);
    std::cout << "LUP factorization: " << std::endl;
    std::cout << "LU: " << std::endl;
    printMatrix(lu.lu, LENGTH, LENGTH);
    std::cout << "P: " << std::endl;
    printMatrix(lu.p, LENGTH, LENGTH);

    double** luSol = solveLU(lu, b);
    std::cout << "LU solution: " << std::endl;
    printVector(luSol);

    return 0;
}
