
#ifndef LINALG_GAUSS_SOLVE_H
#define LINALG_GAUSS_SOLVE_H

#include "matrix.h"
#include <cmath>

const long double G_EPS = 1E-4;

template<typename T>
std::pair<Matrix<T>, Matrix<T>> StraightRun(Matrix<T> A, Matrix<T> B, T& det, std::vector<size_t> &v, size_t & rank) {
    v = {};
    size_t number_of_null_rows = 0;
    for (size_t i0 = 0; i0 < A.HorizDim(); i0++) {
        v.push_back(i0);
    }
    if (B.VertDim() != A.VertDim()) {
        throw BadMatrixDimension();
    }
    for (size_t i = 0; i < min(A.HorizDim(), A.VertDim() - number_of_null_rows); i++) {
        if (std::abs(A(i, i)) < G_EPS) {         ///Перестановка строк
            for (size_t i1 = i; i1 < A.VertDim(); i1++) { //// mozhno vichest - number_of_null_rows
                if (std::abs(A(i1, i)) > G_EPS) {
                    for (size_t j1 = i; j1 < A.HorizDim(); j1++) {
                        std::swap(A(i1, j1), A(i, j1));
                    }
                    for (size_t j2 = 0; j2 < B.HorizDim(); j2++) {
                        std::swap(B(i1, j2), B(i, j2));
                    }
                    det = -det;
                    break;
                }
            }
        }
        if (std::abs(A(i, i)) < G_EPS) { ///нужно переставлять столбцы
            for (size_t j3 = i; j3 < A.HorizDim(); j3++) {
                if (std::abs(A(i, j3)) > G_EPS) {
                    for (size_t i3 = 0; i3 < A.VertDim(); i3++) {
                        std::swap(A(i3, i), A(i3, j3));
                    }
                    std::swap(v[i], v[j3]);
                    break;
                }
            }
        }
        if (std::abs(A(i, i)) > G_EPS)         ///Иначе вся строка нулевая
        {
            T a = A(i, i);
            for (size_t j = i; j < A.HorizDim(); j++) {    ///разделили i-тую строку на A[i, i]
                A(i, j) = A(i, j) / a;
            }
            det = det * a;
            for (size_t j2 = 0; j2 < B.HorizDim(); j2++) {
                B(i, j2) = B(i, j2) / a;
            }
            for (size_t k = i + 1; k < A.VertDim(); k++) { /// из всех остальных строк вычли i-тую на A[k, i]
                T c = A(k, i);
                for (size_t j = i; j < A.HorizDim(); j++) {
                    A(k, j) = A(k, j) - c * A(i, j);
                }
                for (size_t j2 = 0; j2 < B.HorizDim(); j2++) {
                    B(k, j2) = B(k, j2) - c * B(i, j2);
                }
            }
        }
        else {     ///вся эта строка нулевая, переставляем на последнее место, про которое точно не известно, что оно нулевое
            number_of_null_rows ++;
            size_t M = A.VertDim() - number_of_null_rows;
            for (size_t j1 = i; j1 < A.HorizDim(); j1++) {
                std::swap(A(M, j1), A(i, j1));
            }
            for (size_t j2 = 0; j2 < B.HorizDim(); j2++) {
                std::swap(B(M, j2), B(i, j2));
            }
            i--;
        }
    }
    for (size_t i = 0; i < min(A.HorizDim(), A.VertDim()); i++) {
        if (std::abs(A(i, i)) > G_EPS) {
            rank++;
        }
    }
    return {A, B};
}

template<typename T>
std::pair<Matrix<T>, Matrix<T>> ReverseRun(Matrix<T> A, Matrix<T> B, std::vector<size_t> &v, size_t& rank) {
    if (B.VertDim() != A.VertDim()) {
        throw BadMatrixDimension();
    }
    for (size_t i = rank-1; i < size_t(-1); i--) {
        for (size_t k = i - 1; k < size_t(-1); k--) {
            T c = A(k, i);
            for (size_t j = i; j < A.HorizDim(); j++) {
                A(k, j) = A(k, j) - c * A(i, j);
            }
            for (size_t j1 = 0; j1 < B.HorizDim(); j1++) {
                B(k, j1) = B(k, j1) - c * B(i, j1);
            }
        }
    }
    return {A, B};
}



#endif //LINALG_GAUSS_SOLVE_H
