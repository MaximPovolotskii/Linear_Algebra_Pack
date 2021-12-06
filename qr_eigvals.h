//
// Created by Максим on 05.12.2021.
//

#ifndef LINALG_QR_EIGVALS_H
#define LINALG_QR_EIGVALS_H

#include "matrix.h"
#include <algorithm>
#include <cmath>

const double THRESHOLD = 1E-3;

//the algorithm works only for symmetric matrices

template <typename T>
bool Converged (Matrix<T>& mat) {
    for (size_t i = 1; i < mat.VertDim(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (abs(mat(i, j)) > THRESHOLD)
                return false;
        }
    }
    return true;
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> QR_Eigenvalues (Matrix<T>& A) {
    bool converged = false;
    Matrix<T> A_k = A;
    Matrix<T> U_k(A.VertDim(), A.HorizDim(), IDENTITY);
    std::pair<Matrix<T>, Matrix<T>> QR_k;

    while (!coverged) {
        QR_k = QRDecomposition(A_k); //Q_k, R_k
        A_k = QR_k.second * QR_k.first; //A_(k+1) = R_k * Q_k
        U_k = U_k * QR_k.first; //U_(k+1) = U_k * Q_k
        converged = Converged(A_k);
    }

    return {A_k, U_k};
}

#endif //LINALG_QR_EIGVALS_H
