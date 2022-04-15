#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <chrono>
#include <climits>

#include "vector.h"
#include "matrix.h"
#include "lu.h"
#include "ldlt.h"

#define VAR 9

//Stores min, max and sum of norms for the report
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

//Generates symmetric matrix as required in p.1
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

//Exports data for the residual-iterations plot
void getDiagramData(double li, double hi, Matrix m, double param, const std::string& fileName) {
    std::ofstream fout(fileName);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    Vector y = Vector(m.length());
    for(int i = 0; i < m.length(); i++) {
        y(i) = distribution(gen);
    }

    Vector b = m * y;

    Sor sor = m.solveSOR(b, param, 1e-8);
    for(int i = 0; i < sor.step; i++) {
        fout << sor.residual[i] << " ";
    }

}

//Explores how vector b affects residual
void exploreRes(double li, double hi, Matrix m, const std::string& fileName) {
    std::ofstream fout(fileName);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    Vector y = Vector(m.length());
    for (int i = 0; i < m.length(); i++) {
        y(i) = distribution(gen);
    }
    Vector b = m * y;

    std::uniform_real_distribution<double> smallDist(-0.1, 0.1);
    Vector uv = Vector(m.length());
    for (int i = 0; i < m.length(); i++) {
        uv(i) = smallDist(gen);
    }

    for (int i = 0; i < 5; i++) {
        Vector bNew = b;
        for (int j = 0; j < m.length(); j++) {
            uv(j) = smallDist(gen);
        }
        if(i != 0) {
            bNew = b + uv;
        }

        Vector gauss = m.solveGauss(bNew);

        Vector gaussCol = m.solveGaussCol(bNew);

        Vector luSol = LU(m).solve(bNew);

        Vector ldltSol = LDLT(m).solve(bNew);

        Sor sor = m.solveSOR(bNew, (1 - (double) 9 / 40), 1e-8);

        if(i == 0) {
            fout << "---------Default Vector b norms-----------" << std::endl;
        }
        else {
            fout << "---------Diff Norm " << i + 1 << "-----------" << std::endl;
        }
        fout << "Gauss: " << (gauss - y).norm() << std::endl;
        fout << "Gauss Column: " << (gaussCol - y).norm() << std::endl;
        fout << "LUP: " << (luSol - y).norm() << std::endl;
        fout << "LDLT: " << (ldltSol - y).norm() << std::endl;
        fout << "SOR: " << (sor.x - y).norm() << std::endl;
        fout << std::endl;
    }

}

//Generates report as required in p.9
void checkMatrix(double li, double hi, Matrix m, const std::string& fileName) {
    std::ofstream fout(fileName);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(li, hi);

    fout << "Initial Matrix: " << m << std::endl;

    Vector y = Vector(m.length());
    for(int i = 0; i < m.length(); i++) {
        y(i) = distribution(gen);
    }
    fout << "Vector y: ";
    fout << y << std::endl;

    Vector b = m * y;
    fout << "Vector b: " << b << std::endl;

    Matrix inverseM = m.inverse();
    fout << "Inverse Matrix: " << inverseM << std::endl;
    fout << "Condition number: " << m.condition() << std::endl;

    Vector gauss = m.solveGauss(b);
    fout << "Gauss: " << gauss << std::endl;

    Vector gaussCol = m.solveGaussCol(b);
    fout << "Gauss column: " << gaussCol << std::endl;

    LU lup = LU(m);
    fout << "LUP factorization: " << std::endl;
    fout << "LU: " << lup.getLU() << std::endl;
    fout << "P: " << lup.getP() << std::endl;

    Vector luSol = lup.solve(b);
    fout << "LUP solution: " << luSol << std::endl;

    Vector ldltSol = LDLT(m).solve(b);

    Sor sor = m.solveSOR(b, (1 - (double)9 / 40), 1e-8);
    fout << "SOR: " << sor.x << std::endl;

    fout << "-----------------NORM(SOLUTION - PRECISE)------------------" << std::endl;
    fout << "Gauss norm: " << (gauss - y).norm() << std::endl;
    fout << "GaussCol norm: " << (gaussCol - y).norm() << std::endl;
    fout << "LUP norm: " << (luSol - y).norm() << std::endl;
    fout << "LDLT norm: " << (ldltSol - y).norm() << std::endl;
    fout << "SOR norm: " << (sor.x - y).norm() << std::endl;

}

//Generates report as required in p.8
void generateReport(int var, int tries, const std::string& reportFileName, const std::string& maxCondFileName) {
    std::ofstream fout(reportFileName);
    std::ofstream mout(maxCondFileName);

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

    mData inverse, gauss, gaussCol, luFactorize, luSolution, ldltSolution, sor;

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
        Vector ldltSol = LDLT(m).solve(b);
        end = std::chrono::high_resolution_clock::now();
        ldltSolution.sumTime += std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

        cur = (ldltSol - y).norm();
        saveDiffNorm(cur, &ldltSolution);


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

    fout << "----------LDLT FACTORIZATION AND SOLUTION----------" << std::endl;
    fout << "Avg time: " << (double)ldltSolution.sumTime / tries << " ms" << std::endl;
    fout << std::endl;
    fout << "Min norm: " << ldltSolution.min << std::endl;
    fout << "Max norm: " << ldltSolution.max << std::endl;
    fout << "Avg norm: " << ldltSolution.sumVal / tries << std::endl;
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

    generateReport(VAR, 100, "reports/report.txt", "reports/maxMatrix.txt");

    auto** arrA1 = new double*[4];
    arrA1[0] = new double[4]{std::pow(VAR, 2) + 15, VAR - 1, -1, -2};
    arrA1[1] = new double[4]{VAR - 1, -15 - std::pow(VAR, 2), -VAR + 4, -4};
    arrA1[2] = new double[4]{-1, -VAR + 4, std::pow(VAR, 2) + 8, -VAR};
    arrA1[3] = new double[4]{-2, -4, -VAR, std::pow(VAR, 2) + 10};

    Matrix a1 = Matrix(4, arrA1);
    checkMatrix(li, hi, a1, "reports/a1.txt");

    Matrix A2 = Matrix(8);
    A2(0,0) = 1;
    for(int i = 1; i < 8; i++){
        A2(0,i) = i + VAR;
    }
    for(int i = 0, j = 100; i < 4; i++, j *= 10){
        A2(1,i) = j * VAR;
    }
    for(int i = 4, j = -1000; i < 7; i++, j *= 10) {
        A2(1, i) = j * VAR;
    }
    A2(1,7) = 1;
    for(int i = 0; i < 8; i++) {
        A2(2, i) = -i + VAR;
    }
    for(int i = 0, j = 1; i < 5; i++, j *= 10) {
        A2(3, i) = j * VAR - 1000;
    }
    for(int i = 0; i < 3; i++) {
        A2(3, i + 5) = -VAR + i;
    }
    A2(4,0) = VAR;
    for(int i = 1; i < 8; i++) {
        A2(4, i) = 1 - i;
    }
    for(int i = 0, j = 2019, k = 1; i < 8; i++, j++, k = -k) {
        A2(5, i) = k * (VAR - j);
    }
    for(int i = 0, j = 2000, k = 2; i < 5; i++, j += 5, k *= 2) {
        A2(6, i) = k * VAR - j;
    }
    A2(6,5) = 2019 * VAR;
    A2(6,6) = -2020 * VAR;
    A2(6,7) = 2021 * VAR;
    A2(7,0) = 1020 - 2 * VAR;
    A2(7,1) = -2924 + 896 * VAR;
    A2(7,2) = 1212 + 9808 * VAR;
    A2(7,3) = -2736 + 98918 * VAR;
    A2(7,4) = 1404 - 11068 * VAR;
    A2(7,5) = -1523 - 8078 * VAR;
    A2(7,6) = 2625 - 102119 * VAR;
    A2(7,7) = -1327 + 1924 * VAR;

    checkMatrix(li, hi, A2.transpose() * A2, "reports/a2.txt");


    std::ifstream fin("reports/maxMatrix.txt");
    Matrix maxCond = Matrix(256);
    for(int i = 0; i < 256; i++) {
        for(int j = 0; j < 256; j++) {
            fin >> maxCond(i, j);
        }
    }

    exploreRes(li, hi, maxCond, "reports/normsMaxCond.txt");

    getDiagramData(li, hi, maxCond, 0.8, "reports/plot/a1Plot08.txt");
    getDiagramData(li, hi, maxCond, 1.0, "reports/plot/a1Plot10.txt");
    getDiagramData(li, hi, maxCond, 1.2, "reports/plot/a1Plot12.txt");


    return 0;
}
