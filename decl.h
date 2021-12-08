#include <stdexcept>

#pragma once

class BadMatrixDimension: public std::exception {
public:
    const char* what() const throw() override {
        return "Dimensions do not converge";
    }
};

class NotSquareMatrix: public BadMatrixDimension {
public:
    const char* what() const throw() override {
        return "Matrix is not square";
    }
};

class IndexOutOfMatrix: public std::exception {
    const char* what() const throw() override {
        return "Index out of range";
    }
};

template<class T>
const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}