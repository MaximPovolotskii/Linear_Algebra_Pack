
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
std::pair<Matrix<T>, Matrix<T>> StraightRun(Matrix<T> A, Matrix<T> B, T& det) {
    if (A.HorizDim() != A.VertDim()) {
        throw BadMatrixDimension(); // пока так
    }
    if (B.VertDim() != A.VertDim()) {
        throw BadMatrixDimension();
    }
    for (size_t i = 0; i < min(A.HorizDim(), A.VertDim()); i++) {
        if (abs(A(i, i)) < EPS) {
            for (size_t i1 = i; i1 < A.VertDim(); i1++) {
                if (A(i1, i) > EPS) {
                    for (size_t j1 = i; j1 < A.HorizDim(); j1++) {
                        std::swap(A(i1, j1), A(i, j1));
                    }
                    for (size_t j2 = 0; j2 < A.HorizDim(); j2++) {
                        std::swap(B(i1, j2), B(i, j2));
                    }
                    det = det * (-1);
                    break;
                }
            }
        }
        if (abs(A(i, i)) > EPS)         ///Проверка на ноль
        {
            T a = A(i, i);
            for (size_t j = i; j < A.HorizDim(); j++) {    ///разделили i-тую строку на A[i, i]
                A(i, j) = A(i, j) / a;
            }
            det = det * a;
            for (size_t j2 = 0; j2 < A.HorizDim(); j2++) {
                B(i, j2) = B(i, j2) / a;
            }
            for (size_t k = i + 1; k < A.VertDim(); k++) { /// из всех остальных строк вычли i-тую на A[k, i]
                T c = A(k, i);
                for (size_t j = i; j < A.HorizDim(); j++) {
                    A(k, j) = A(k, j) - c * A(i, j);
                }
                for (size_t j2 = 0; j2 < A.HorizDim(); j2++) {
                    B(k, j2) = B(k, j2) - c*B(i, j2);
                }
            }
        }
    }
    return {A, B};
}

template<typename T>
std::pair<Matrix<T>, Matrix<T>> ReverseRun(Matrix<T> A, Matrix<T> B) {
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
                for (size_t j1 = 0; j1 < A.HorizDim(); j1++) {
                    B(k, j1) = B(k, j1) - c * B(i, j1);
                }
            }
        }
    }
    return {A, B};
}



#endif //LINALG_GAUSS_SOLVE_H
