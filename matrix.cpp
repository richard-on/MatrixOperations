#include <iostream>
#include <iomanip>
#include <climits>

#include "matrix.h"

Matrix::Matrix(int len) {
    if(len < 2) {
        throw std::invalid_argument("incorrect Matrix dimensions");
    }

    this->len = len;
    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len]{0};
    }
}

Matrix::Matrix(int len, double **data) {
    if(len < 2) {
        throw std::invalid_argument("incorrect Matrix dimensions");
    }

    this->len = len;
    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len];
        for(int j = 0; j < len; j++) {
            this->data[i][j] = data[i][j];
        }
    }
}

Matrix::Matrix(const Matrix &other) {
    this->len = other.len;

    this->data = new double*[len];
    for(int i = 0; i < len; i++) {
        this->data[i] = new double[len];
        for(int j = 0; j < len; j++) {
            this->data[i][j] = other(i, j);
        }
    }
}

Matrix::Matrix(Matrix &&other) noexcept
    : data(other.data), len(other.len)
    {

    other.data = nullptr;
    other.len = 0;
}


Matrix &Matrix::operator=(const Matrix& other) {
    if(this != &other){
        for(int i = 0; i < len; i++) {
            delete[] data[i];
        }
        delete[] data;

        this->len = other.len;
        this->data = new double*[len];
        for(int i = 0; i < len; i++) {
            this->data[i] = new double[len];
            for(int j = 0; j < len; j++) {
                this->data[i][j] = other.data[i][j];
            }
        }
    }

    return *this;
}

Matrix &Matrix::operator=(Matrix &&other) noexcept {
    if (this != &other) {
        for (int i = 0; i < len; i++) {
            delete[] data[i];
        }
        delete[] data;

        this->len = other.len;
        this->data = new double*[len];
        for (int i = 0; i < len; i++) {
            this->data[i] = new double[len];
            for (int j = 0; j < len; j++) {
                this->data[i][j] = other.data[i][j];
            }
        }

        other.data = nullptr;
        other.len = 0;
    }

    return *this;
}

double& Matrix::operator()(int row, int col) const {
    return this->data[row][col];
}

bool Matrix::operator==(const Matrix &other) const {
    if(this->len == other.len){
        for(int i = 0; i < len; i++) {
            for(int j = 0; j < len; j++) {
                if(this->data[i][j] != other.data[i][j]) {
                    return false;
                }
            }
        }
        return true;
    }
    return false;
}

bool Matrix::operator!=(const Matrix &other) const {
    if(this == &other){
        return false;
    }
    return true;
}

Matrix Matrix::operator+(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            m.data[i][j] = this->data[i][j] + other.data[i][j];
        }
    }

    return m;
}

Matrix Matrix::operator-(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            m.data[i][j] = this->data[i][j] - other.data[i][j];
        }
    }

    return m;
}

Matrix Matrix::operator*(const Matrix &other) const {
    if(this->len != other.len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Matrix m = Matrix(len);
    for(int i = 0; i < len; i++) {
        for(int j = 0; j < len; j++) {
            for(int k = 0; k < len; k++) {
                m.data[i][j] += this->data[i][k] * other.data[k][j];
            }
        }
    }

    return m;
}

Vector Matrix::operator*(const Vector &other) const {
    if(other.length() != this->len) {
        throw std::invalid_argument("error trying to compare Matrix of different length");
    }

    Vector v = Vector(this->len);
    for(int i = 0; i < v.length(); i++) {
        for(int k = 0; k < v.length(); k++) {
            v(i) += this->data[i][k] * other(k);
        }
    }

    return v;
}

std::ostream &operator<<(std::ostream &ostream, const Matrix &matrix) {
    ostream << std::endl;
    for (int i = 0; i < matrix.length(); i++) {
        if(i == 0) {
            ostream << "/";
        }
        else if(i == matrix.length() - 1) {
            ostream << "\\";
        }
        else{
            ostream << "|";
        }

        for (int j = 0; j < matrix.length(); j++) {

            ostream << std::setw(8) << std::setprecision(3) << matrix(i, j) << " ";
        }
        if(i == 0) {
            ostream << "\\";
        }
        else if(i == matrix.length() - 1) {
            ostream << "/";
        }
        else{
            ostream << "|";
        }
        ostream << std::endl;

    }

    return ostream;
}


double **Matrix::getData() const {
    return data;
}

int Matrix::length() const {
    return len;
}

double Matrix::norm() {
    double norm = INT_MIN;
    double columnSum = 0;

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            columnSum += std::abs(this->data[i][j]);
        }
        if(columnSum > norm) {
            norm = columnSum;
        }
        columnSum = 0;
    }

    return norm;
}

double Matrix::condition() {
    double norm = this->norm();
    double invNorm = this->inverse().norm();

    return norm * invNorm;
}

Matrix Matrix::inverse() {
    auto res = Matrix(len);
    auto** m = new double*[len];
    for (int i = 0; i < len; i++) {
        m[i] = new double[len * 2]{0};
        for (int j = 0; j < 2 * len; j++) {
            if(j < len) {
                m[i][j] = this->data[i][j];
            }
            if (j == (i + len))
                m[i][j] = 1;
        }
    }

    for (int i = len-1; i > 0; i--) {
        if (m[i - 1][0] < m[i][0]) {
            double* temp = m[i];
            m[i] = m[i - 1];
            m[i - 1] = temp;
        }
    }

    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            if (j != i) {
                double d = m[j][i] / m[i][i];
                for (int k = 0; k < 2 * len; k++) {
                    m[j][k] -= m[i][k] * d;
                }
            }
        }
    }

    for (int i = 0; i < len; i++) {
        double d = m[i][i];
        for (int j = 0; j < 2 * len; j++) {
            m[i][j] = m[i][j] / d;
        }
    }

    for (int i = 0; i < len; i++) {
        int k = 0;
        for (int j = len; j < 2 * len; j++) {
            res(i, k) = m[i][j];
            k++;
        }
    }

    for (int i = 0; i < len; i++) {
        delete[] m[i];
    }
    delete[] m;

    return res;
}


Matrix::~Matrix() {
    for (int i = 0; i < len; i++) {
        delete[] data[i];
    }
    delete[] data;
}
