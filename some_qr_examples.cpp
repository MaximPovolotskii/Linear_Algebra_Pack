#include <iostream>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "gauss_run.h"
#include "qr_decomposition.h"
#include "qr_eigvals.h"
#include <tuple>
#include <complex>

using namespace std::complex_literals;
int main() {

    std::complex<long double> a[16] {0.f+1.if, -1, 1, 0, 0, 6, 3, 3, -4, 8.if, 2.f+1.if, 5, -1, 0, 6, 0};

    Matrix<std::complex<long double>> A(a, 4, 4);
    Matrix<std::complex<long double>> Q(a, 4, 4);
    Matrix<std::complex<long double>> R(a, 4, 4);

    std::pair<std::vector<std::complex<long double>>, Matrix<long double>> pqr(QREigenvalues(A));
    std::vector<std::complex<long double>> im = pqr.first;
    Matrix<long double> M = pqr.second;

    std::cout << "eigenvals"<<"\n";
    for (std::complex<long double> i : im) {
        std::cout << i <<" ";
    }
    std::cout << '\n';


    std::cout << "eigenvectors";
    std::cout << '\n';
    for (size_t i = 0; i < M.VertDim(); i++) {
        for (size_t j = 0; j < M.HorizDim(); j++) {
            std::cout << M(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';*/


    QRDecomposition(A, Q, R);
    std::cout << 'Q';
    std::cout << '\n';
    for (size_t i = 0; i < Q.VertDim(); i++) {
        for (size_t j = 0; j < Q.HorizDim(); j++) {
            std::cout << Q(i, j) << " ";
        }
    std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';

    std::cout << 'R';
    std::cout << '\n';
    for (size_t i = 0; i < R.VertDim(); i++) {
        for (size_t j = 0; j < R.HorizDim(); j++) {
            std::cout << R(i, j) << " ";
        }
        std::cout << '\n';
    }



    p = QREigenvalues(A);

   std::cout<<"M1"<<"\n";
    Matrix<long double> M1 = p.first;
    for (size_t i = 0; i < M1.VertDim(); i++) {
        for (size_t j = 0; j < M1.HorizDim(); j++) {
            std::cout << M1(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';

    Matrix<long double> M2 = p.second;
    for (size_t i = 0; i < M2.VertDim(); i++) {
        for (size_t j = 0; j < M2.HorizDim(); j++) {
            std::cout << M2(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';


    auto Q1 = Q;
    Q1.Conjugate();
    Q1.Transpose();


}
