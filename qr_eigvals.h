
#ifndef LINALG_QR_EIGVALS_H
#define LINALG_QR_EIGVALS_H

#include "matrix.h"
#include <algorithm>
#include <cmath>
#include "qr_decomposition.h"

const double THRESHOLD = 1E-6;
const size_t MAXITER = 1E8;

template <typename T>
bool Converged (const Matrix<T>& mat, double threshold) {
    for (size_t i = 2; i < mat.VertDim(); ++i) {
        for (size_t j = 0; j < i - 1; ++j) {
            if (std::abs(mat(i, j)) > threshold)
                return false;
        }
    }
    for (size_t i = 1; i < mat.VertDim() - 1; ++i) {
        if (std::abs(mat(i, i-1)) > threshold && std::abs(mat(i + 1, i)) > threshold)
            return false;
    }
    return true;
}


template <typename T>
std::pair<std::complex<T>, std::complex<T>> ComputeComplexEigenvalues(const T& a11, const T& a12,
                                                                      const T& a21, const T& a22) {
    T D = (a11 + a22)*(a11 + a22) - 4 * (a11 * a22 - a12 * a21);
    if (D < 0) {
        T re = (a11 + a22)/2;
        T im = std::sqrt(a11 * a22 - a12 * a21 - re*re);
        return {std::complex<T>(re, im), std::complex<T>(re, -im)};
    }
    else {
        T x1 = (a11 + a22 + std::sqrt(D))/2;
        T x2 = (a11 + a22 - std::sqrt(D))/2;
        return {std::complex<T>(x1, 0), std::complex<T>(x2, 0)};
    }
}


template <typename T>
std::pair<std::vector<std::complex<T>>, Matrix<T>> QREigenvalues (const Matrix<T>& A,
                                                                                   size_t max_iterations = MAXITER,
                                                                                   double threshold = THRESHOLD) {
    if (A.VertDim() != A.HorizDim()) {
        std::cerr<<"In QREigenvalues (const Matrix<T>& )    ";
        throw NotSquareMatrix();
    }
    bool converged = false;
    Matrix<T> A_k = A;
    Matrix<T> U_k(A.VertDim(), A.HorizDim(), IDENTITY);
    Matrix<T> Q_k(A.VertDim(), A.VertDim());
    Matrix<T> R_k(A.VertDim(), A.HorizDim());
    size_t count = 0;
    while (!converged && count < max_iterations) {
        QRDecomposition(A_k, Q_k, R_k); //Q_k, R_k

        A_k = R_k * Q_k; //A_(k+1) = R_k * Q_k
        U_k = U_k * Q_k; //U_(k+1) = U_k * Q_k
        converged = Converged(A_k, threshold);
        ++count;
    }
    std::vector<std::complex<T>> img;
    std::pair<std::complex<T>,  std::complex<T>> p_img;


    for (size_t i = 0; i < A_k.VertDim()-1; i++) {
        if (std::abs(A_k(i+1, i)) > threshold) {
            p_img = ComputeComplexEigenvalues(A_k(i, i), A_k(i, i+1), A_k(i+1, i), A_k(i+1, i+1));
            img.push_back(p_img.first);
            img.push_back(p_img.second);
            i++;
        }
        else {
            img.push_back(std::complex<T>(A_k(i, i), 0));
        }
    }
    if (std::abs(A_k(A_k.VertDim()-1, A_k.VertDim()-2)) < threshold) {
        img.push_back(A_k(A_k.VertDim()-1, A_k.VertDim()-1));
    }
    return {img, U_k};
}

#endif //LINALG_QR_EIGVALS_H
