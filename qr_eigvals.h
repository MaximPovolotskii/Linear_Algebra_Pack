
#ifndef LINALG_QR_EIGVALS_H
#define LINALG_QR_EIGVALS_H

#include "matrix.h"
#include <algorithm>
#include <cmath>
#include "qr_decomposition.h"
#include "linear_system_solve.h"

const double THRESHOLD = 1E-8;
const size_t MAXITER = 1E10;

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
std::pair<std::complex<T>, std::complex<T>> ComputeComplexEigenvalues2x2(const T& a11, const T& a12,
                                                                      const T& a21, const T& a22, bool& are_complex) {
    T D = (a11 + a22)*(a11 + a22) - 4 * (a11 * a22 - a12 * a21);
    if (D < 0) {
        T re = (a11 + a22)/2;
        T im = std::sqrt(a11 * a22 - a12 * a21 - re*re);
        are_complex = true;
        return {std::complex<T>(re, im), std::complex<T>(re, -im)};
    }
    else {
        T x1 = (a11 + a22 + std::sqrt(D))/2;
        T x2 = (a11 + a22 - std::sqrt(D))/2;
        are_complex = false;
        return {std::complex<T>(x1, 0), std::complex<T>(x2, 0)};
    }
}


template <typename T>
std::pair<std::vector<T>, std::vector<std::complex<T>>> QREigenvalues (const Matrix<T>& A,
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
    std::vector<T> real;
    std::vector<std::complex<T>> img;
    std::pair<std::complex<T>,  std::complex<T>> p_img;

    for (size_t i = 0; i < A_k.VertDim()-1; i++) {
        if (std::abs(A_k(i+1, i)) > threshold) {
            bool are_complex;
            p_img = ComputeComplexEigenvalues2x2(A_k(i, i), A_k(i, i+1), A_k(i+1, i), A_k(i+1, i+1), are_complex);
            if (are_complex) {
                img.push_back(p_img.first);
                img.push_back(p_img.second);}
            else{
                real.push_back(std::real(p_img.first));
                real.push_back(std::real(p_img.second));
            }
            i++;
        }
        else {
            real.push_back(A_k(i, i));
        }
    }
    if (std::abs(A_k(A_k.VertDim()-1, A_k.VertDim()-2)) < threshold) {
        real.push_back(A_k(A_k.VertDim()-1, A_k.VertDim()-1));
    }
    return {real, img};
}

template <typename T>
std::pair<size_t, Matrix<T>> FindEigenvectors(const Matrix<T>& A, T lambda) {
    if (A.VertDim() != A.HorizDim()) {
        std::cerr<<"In FindEigenvector (const Matrix<T>& , T)    ";
        throw NotSquareMatrix();
    }
    Matrix<T> B = (A - Matrix<T>(A.HorizDim(), A.VertDim(), IDENTITY)*lambda);

    std::tuple<bool, Vector<T>, Matrix<T>> solution = LinearSystemSolve(B, Vector<T>(B.VertDim(), NULLMATRIX));
    if (!std::get<0>(solution)) {
        return {0, Vector<T>(0)};
    }
    else {

        return {std::get<2>(solution).HorizDim(), std::get<2>(solution)};
    }
}

template <typename T>
std::pair<size_t, Matrix<std::complex<T>>> FindEigenvectors(const Matrix<T>& A, std::complex<T> lambda) {
    Matrix<std::complex<T>> A_c(A.VertDim(), A.HorizDim());
    for (size_t i = 0; i < A_c.VertDim(); i++) {
        for (size_t j = 0; j < A_c.HorizDim(); j++) {
            A_c(i, j) = std::complex<T>(A(i, j));
        }
    }
    return FindEigenvectors(A_c, lambda);
}

/*template <typename T>
std::vector<size_t, Matrix<std::complex<T>>> FindAllEigenvectors(const Matrix<T>& A) {

}*/

template <typename T>
std::tuple<size_t, Vector<T>, Matrix<T>> FindGeneralizedEigenvectors(const Matrix<T>& A, std::complex<T> lambda,
                                                                     const Vector<T>& v ) {
    if (A.VertDim() != A.HorizDim()) {
        std::cerr<<"In FindGeneralizedEigenvectors (const Matrix<T>& , T)    ";
        throw NotSquareMatrix();
    }

    Matrix<T> B = (A - Matrix<T>(A.HorizDim(), A.VertDim(), IDENTITY)*lambda);

    std::tuple<bool, Vector<T>, Matrix<T>> solution = LinearSystemSolve(B, v);
    if (!std::get<0>(solution)) {
        return {0, Vector<T>(0), Vector<T>(0, 0)};
    }
    else {
        return {std::get<2>(solution).HorizDim(),std::get<1>(solution), std::get<2>(solution)};
    }
}


#endif //LINALG_QR_EIGVALS_H
