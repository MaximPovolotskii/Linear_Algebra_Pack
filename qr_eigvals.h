
#ifndef LINALG_QR_EIGVALS_H
#define LINALG_QR_EIGVALS_H

#include "matrix.h"
#include <algorithm>
#include <cmath>
#include "qr_decomposition.h"

const double THRESHOLD = 1E-1;


template <typename T>
bool Converged (Matrix<T>& mat) {
    for (size_t i = 1; i < mat.VertDim(); ++i) {
        for (size_t j = 0; j < i; ++j) {
            if (std::abs(mat(i, j)) > THRESHOLD)
                return false;
        }
    }
    return true;
}

template <typename T>
std::pair<Matrix<T>, Matrix<T>> QREigenvalues (const Matrix<T>& A) {
    bool converged = false;

    Matrix<T> A_k = A;
    Matrix<T> U_k(A.VertDim(), A.HorizDim(), IDENTITY);
    Matrix<T> Q_k(A.VertDim(), A.VertDim());
    Matrix<T> R_k(A.VertDim(), A.HorizDim());
    size_t i = 0;

    while (i < 1200) {
        QRDecomposition(A_k, Q_k, R_k); //Q_k, R_k
        A_k = R_k * Q_k; //A_(k+1) = R_k * Q_k
        U_k = U_k * Q_k; //U_(k+1) = U_k * Q_k
        converged = Converged(A_k);
        ++i;
    }

    return {A_k, U_k};
}

#endif //LINALG_QR_EIGVALS_H
