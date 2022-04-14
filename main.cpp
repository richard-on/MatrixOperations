#include <iostream>
#include <fstream>
#include <random>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <climits>

#include "vector.h"
#include "matrix.h"
#include "lu.h"

#define VAR 9
#define LENGTH 3

/*void report() {
    std::ofstream fout("report.txt");

    double high_interval = std::pow(2, double(VAR)/4);
    double low_interval = -high_interval;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> distribution(low_interval, high_interval);

    double minCond = INT_MAX;
    double maxCond = INT_MIN;
    double** maxCondMatrix;
    double avgCond = 0;

    double avgInverseTime = 0;

    for(int t = 0; t < 100; t++) {
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

        double** b = multiplyMatrix(matrix, y, LENGTH, 1);

        std::chrono::high_resolution_clock::time_point begin = std::chrono::high_resolution_clock::now();
        double** inverse = inverseMatrix(matrix);
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        avgInverseTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();

        double mc = cubicNorm(matrix);
        double mi = cubicNorm(inverse);
        double cond = mc * mi;
        avgCond += cond;
        if(cond < minCond) {
            minCond = cond;
        }
        if(cond > maxCond) {
            maxCond = cond;
            maxCondMatrix = matrix;
        }

        //double** gauss = solveGauss(Matrix, b);

        double** gaussCol = solveGaussCol(matrix, b);

        auto** newMatrix = new double*[LENGTH];
        for(int i = 0; i < LENGTH; i++) {
            newMatrix[i] = new double[LENGTH];
            for(int j = 0; j < LENGTH; j++) {
                newMatrix[i][j] = matrix[i][j];
            }
        }
        LU lu = luFactorize(matrix);

        double** luSol = solveLU(lu, b);

        double** sor = sorSolve(newMatrix, b, (1 - (double)9/40));
    }

    fout << "Average condition: " << avgCond / 100 << std::endl;
    fout << "Min condition: " << minCond << std::endl;
    fout << "Max condition: " << maxCond << std::endl;
    fout << "Max condition Matrix: " << maxCondMatrix[0][0] << std::endl;

    fout << "Average inverse time: " << avgInverseTime / 100 << " ns" << std::endl;
}*/

Matrix generateMatrix(int length, double li, double hi){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    double rowSum = 0;
    Matrix m = Matrix(length);
    for(int i = 0; i < m.length(); i++) {
        for(int j = m.length() - 1; j >= 0; j--) {
            if(j > i) {
                m(i, j) = distribution(gen);
            }
            else if(j < i) {
                m(i, j) = m(j, i);
            }
            rowSum += std::abs(m(i, j));

            if(j == 0){
                m(i, i) = rowSum + 1;
                rowSum = 0;
            }
        }
    }

    return m;
}

void example(double li, double hi) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> distribution((int(li)), int(hi));

    Matrix m = generateMatrix(5, li, hi);
    std::cout << "Initial Matrix: " << m << std::endl;

    Vector y = Vector(LENGTH);
    for(int i = 0; i < m.length(); i++) {
        y(i) = distribution(gen);
    }
    std::cout << "Vector y: " << y << std::endl;

    Vector b = m * y;
    std::cout << "Vector b: " << b << std::endl;

    Matrix inverseM = m.inverse();
    std::cout << "Inverse Matrix: " << inverseM << std::endl;
    std::cout << "Condition number: " << m.condition() << std::endl;

    Vector gauss = m.solveGauss(b);
    std::cout << "Gauss: " << gauss << std::endl;

    Vector gaussCol = m.solveGaussCol(b);
    std::cout << "Gauss column: " << gaussCol << std::endl;

    LU lup = LU(m);
    std::cout << "LUP factorization: " << std::endl;
    std::cout << "LU: " << lup.getLU() << std::endl;
    std::cout << "P: " << lup.getP() << std::endl;

    Vector luSol = lup.solve(b);
    std::cout << "LUP solution: " << luSol << std::endl;

    Sor sor = m.solveSOR(b, (1 - (double)9 / 40));
    std::cout << "SOR: " << sor.x << std::endl;
}

struct mData {
    double max = INT_MIN;
    double min = INT_MAX;
    double sumVal = 0;
    long long sumTime = 0;
};

void saveDiffNorm(double cur, mData* method) {
    if(cur < method->min) {
        method->min = cur;
    }
    if(cur > method->max) {
        method->max = cur;
    }
    method->sumVal += cur;
}

void generateReport(int var, int tries) {
    std::ofstream fout("report.txt");
    std::ofstream mout("maxMatrix.txt");

    double hi = std::pow(2, double(VAR)/4);
    double li = -hi;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    double minCond = INT_MAX;
    double maxCond = INT_MIN;
    double sumCond = 0;
    Matrix maxCondMatrix;

    double minStep = INT_MAX;
    double maxStep = INT_MIN;
    double sumStep = 0;

    mData inverse, gauss, gaussCol, luFactorize, luSolution, sor;

    for(int t = 0; t < tries; t++) {
        Matrix m = generateMatrix(256, li, hi);

        Vector y = Vector(256);
        for(int i = 0; i < m.length(); i++) {
            y(i) = distribution(gen);
        }

        Vector b = m * y;

        double cur = m.condition();
        if(cur < minCond) {
            minCond = cur;
        }
        if(cur > maxCond) {
            maxCond = cur;
            maxCondMatrix = m;
        }
        sumCond += cur;

        auto start = std::chrono::high_resolution_clock::now();
        Matrix inverseM = m.inverse();
        auto end = std::chrono::high_resolution_clock::now();
        inverse.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


        start = std::chrono::high_resolution_clock::now();
        Vector gaussSol = m.solveGauss(b);
        end = std::chrono::high_resolution_clock::now();
        gauss.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        cur = (gaussSol - y).norm();
        saveDiffNorm(cur, &gauss);


        start = std::chrono::high_resolution_clock::now();
        Vector gaussColSol = m.solveGaussCol(b);
        end = std::chrono::high_resolution_clock::now();
        gaussCol.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        cur = (gaussColSol - y).norm();
        saveDiffNorm(cur, &gaussCol);


        start = std::chrono::high_resolution_clock::now();
        LU lu = LU(m);
        end = std::chrono::high_resolution_clock::now();
        luFactorize.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();


        start = std::chrono::high_resolution_clock::now();
        Vector luSol = lu.solve(b);
        end = std::chrono::high_resolution_clock::now();
        luSolution.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        cur = (luSol - y).norm();
        saveDiffNorm(cur, &luSolution);


        start = std::chrono::high_resolution_clock::now();
        Sor sorSol = m.solveSOR(b, (1 - (double)var / 40), 1e-8);
        end = std::chrono::high_resolution_clock::now();
        sor.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        cur = (sorSol.x - y).norm();
        saveDiffNorm(cur, &sor);
        if(sorSol.step < minStep) {
            minStep = sorSol.step;
        }
        if(sorSol.step > maxStep) {
            maxStep = sorSol.step;
        }
        sumStep += sorSol.step;
    }
    fout << "----------CONDITIONS----------" << std::endl;
    fout << "Avg condition number: " << (double)sumCond / tries << std::endl;
    fout << "Min condition number: " << minCond << std::endl;
    fout << "Max condition number: " << maxCond << std::endl;
    fout << std::endl;

    for(int i = 0; i < maxCondMatrix.length(); i++) {
        for (int j = 0; j < maxCondMatrix.length(); j++) {
            mout << maxCondMatrix(i,j) << " ";
        }
        mout << std::endl;
    }


    fout << "----------INVERSE----------" << std::endl;
    fout << "Avg inverse time: " << (double)inverse.sumTime / tries << " ms" << std::endl;
    fout << std::endl;

    fout << "----------GAUSS----------" << std::endl;
    fout << "Avg time: " << (double)gauss.sumTime / tries << " ms" << std::endl;
    fout << std::endl;
    fout << "Min norm: " << gauss.min << std::endl;
    fout << "Max norm: " << gauss.max << std::endl;
    fout << "Avg norm: " << gauss.sumVal / tries << std::endl;
    fout << std::endl;

    fout << "----------GAUSS COLUMN----------" << std::endl;
    fout << "Avg time: " << (double)gaussCol.sumTime / tries << " ms" << std::endl;
    fout << std::endl;
    fout << "Min norm: " << gaussCol.min << std::endl;
    fout << "Max norm: " << gaussCol.max << std::endl;
    fout << "Avg norm: " << gaussCol.sumVal / tries << std::endl;
    fout << std::endl;

    fout << "----------LUP FACTORIZATION----------" << std::endl;
    fout << "Avg time: " << (double)luFactorize.sumTime / tries << " ms" << std::endl;
    fout << std::endl;

    fout << "----------LUP SOLUTION----------" << std::endl;
    fout << "Avg time: " << (double)luSolution.sumTime / tries << " ms" << std::endl;
    fout << std::endl;
    fout << "Min norm: " << luSolution.min << std::endl;
    fout << "Max norm: " << luSolution.max << std::endl;
    fout << "Avg norm: " << luSolution.sumVal / tries << std::endl;
    fout << std::endl;

    fout << "----------SOR----------" << std::endl;
    fout << "Avg time: " << (double)sor.sumTime / tries << " ms"  << std::endl;
    fout << std::endl;
    fout << "Min norm: " << sor.min << std::endl;
    fout << "Max norm: " << sor.max << std::endl;
    fout << "Avg norm: " << sor.sumVal / tries << std::endl;
    fout << std::endl;
    fout << "Min steps: " << minStep << std::endl;
    fout << "Max steps: " << maxStep << std::endl;
    fout << "Avg steps: " << (double)sumStep / tries << std::endl;
    fout << std::endl;

}

int main() {
    double hi = std::pow(2, double(VAR)/4);
    double li = -hi;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    generateReport(VAR, 10);

    double arrA1[4][4] = {{std::pow(VAR, 2) + 15, VAR - 1, -1, -2},
                          {VAR - 1, -15 - std::pow(VAR, 2), -VAR + 4, -4},
                          {-1, -VAR + 4, std::pow(VAR, 2) + 8, -VAR},
                          {-2, -4, -VAR, std::pow(VAR, 2) + 10}};

    Matrix a1 = Matrix(4, arrA1);

    std::cout << "Initial Matrix: " << a1 << std::endl;

    Vector y = Vector(4);
    for(int i = 0; i < a1.length(); i++) {
        y(i) = distribution(gen);
    }
    std::cout << "Vector y: " << y << std::endl;

    Vector b = a1 * y;
    std::cout << "Vector b: " << b << std::endl;

    Matrix inverseM = a1.inverse();
    std::cout << "Inverse Matrix: " << inverseM << std::endl;
    std::cout << "Condition number: " << a1.condition() << std::endl;

    Vector gauss = a1.solveGauss(b);
    std::cout << "Gauss: " << gauss << std::endl;

    Vector gaussCol = a1.solveGaussCol(b);
    std::cout << "Gauss column: " << gaussCol << std::endl;

    LU lup = LU(a1);
    std::cout << "LUP factorization: " << std::endl;
    std::cout << "LU: " << lup.getLU() << std::endl;
    std::cout << "P: " << lup.getP() << std::endl;

    Vector luSol = lup.solve(b);
    std::cout << "LUP solution: " << luSol << std::endl;

    Sor sor = a1.solveSOR(b, (1 - (double)9 / 40));
    std::cout << "SOR: " << sor.x << std::endl;




    return 0;
}
