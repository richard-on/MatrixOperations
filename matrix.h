#ifndef MV_MATRIX_H
#define MV_MATRIX_H

#include <iosfwd>
#include "vector.h"

class Matrix {
public:
    explicit Matrix(int len = 2);

    Matrix(int len, double** data);

    Matrix(const Matrix& other);

    Matrix(Matrix&& other) noexcept;


    Matrix& operator = (const Matrix& other);

    Matrix& operator = (Matrix &&other) noexcept;

    double& operator () (int rowIdx, int colIdx) const;

    bool operator == (const Matrix &other) const;

    bool operator != (const Matrix &other) const;

    Matrix operator + (const Matrix &other) const;

    Matrix operator - (const Matrix &other) const;

    Matrix operator * (const Matrix &other) const;

    Vector operator * (const Vector &other) const;

    friend std::ostream& operator << (std::ostream& ostream, const Matrix& matrix);


    double **getData() const;

    int length() const;

    Matrix inverse();

    double norm();

    double condition();


    virtual ~Matrix();

private:
    double** data{};
    int len{};
};


#endif //MV_MATRIX_H
