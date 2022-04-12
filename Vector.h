#ifndef MV_VECTOR_H
#define MV_VECTOR_H

class Vector {
public:
    explicit Vector(int len);
    Vector(const Vector& other);
    double* operator - (const Vector& a);
    void print();
    virtual ~Vector();

private:
    double* data{};
    int len{};
};


#endif //MV_VECTOR_H
