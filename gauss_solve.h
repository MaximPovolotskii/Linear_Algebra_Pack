
#ifndef LINALG_GAUSS_SOLVE_H
#define LINALG_GAUSS_SOLVE_H

const double EPS = 1E-9;

template<class T>
const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}

#include "matrix.h"
template<typename T>
std::pair<Matrix<T>, Vector<T>> StraightRun(const Matrix<T> &A0, const Vector<T> & b0, T& det) {
    Matrix<T> A = A0;
    Vector<T> b = b0;
    if (A.HorizDim() != A.VertDim()) {
        throw BadMatrixDimension(); // пока так
    }
    for (size_t i = 0; i < min(A.HorizDim(), A.VertDim()); i++) {
        if (abs(A(i, i)) < EPS) {
            for (size_t i1 = i; i1 < A.VertDim(); i1++) {
                if (A(i1, i) > EPS) {
                    for (size_t j1 = i; j1 < A.HorizDim(); j1++) {
                        std::swap(A(i1, j1), A(i, j1));
                    }
                    det = det * (-1);
                    std::swap(b[i], b[i1]);
                    break;
                }
            }
        }
        if (abs(A(i, i)) > EPS)
        {
            T a = A(i, i);
            for (size_t j = i; j < A.HorizDim(); j++) {    ///разделили i-тую строку на A[i, i]
                A(i, j) = A(i, j) / a;
            }
            det = det * a;
            b[i] = b[i] / a;
            for (size_t k = i + 1; k < A.VertDim(); k++) { /// из всех остальных строк вычли i-тую на A[k, i]
                T c = A(k, i);
                for (size_t j = i; j < A.HorizDim(); j++) {
                    A(k, j) = A(k, j) - c * A(i, j);
                }
                b[k] = b[k] - c * b[i];
            }
        }
    }
    return {A, b};
}

template<typename T>
std::pair<Matrix<T>, Vector<T>> ReverseRun(const Matrix<T> &A0, const Vector<T> & b0) {
    Matrix<T> A = A0;
    Vector<T> b = b0;
    if (A.HorizDim() != A.VertDim()) {
        throw BadMatrixDimension(); // пока так
    }
    for (size_t i = min(A.HorizDim(), A.VertDim()) - 1; i < size_t(-1); i--) {
        if (A(i, i) != 0) {
            for (size_t k = i - 1; k < size_t(-1); k--) {
                T c = A(k, i);
                for (size_t j = i; j < A.HorizDim(); j++) {
                    A(k, j) = A(k, j) - c * A(i, j);
                }
                b[k] = b[k] - c * b[i];
            }
        }
    }
    return {A, b};
}


#endif //LINALG_GAUSS_SOLVE_H
