//
// Created by Richard on 12.04.2022.
//

#include <iostream>
#include <iomanip>
#include "Vector.h"

Vector::Vector(int len) {
    data = new double[len]{0};
}

Vector::Vector(const Vector &other) {

}


double *Vector::operator-(const Vector& a) {
    return nullptr;
}

void Vector::print() {
    for (int i = 0; i < len; i++) {
        if(i == 0) {
            std::cout << "/";
        }
        else if(i == len - 1) {
            std::cout << "\\";
        }
        else{
            std::cout << "|";
        }
        std::cout << std::setw(6) << std::setprecision(3) << data[i];
        if(i == 0) {
            std::cout << "\\";
        }
        else if(i == len - 1) {
            std::cout << "/";
        }
        else{
            std::cout << "|";
        }

        std::cout << std::endl;
    }

    std::cout << std::endl;
}

Vector::~Vector() {
    delete[] data;
}




