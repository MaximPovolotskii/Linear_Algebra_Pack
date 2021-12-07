#include <iostream>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "gauss_solve.h"
#include "qr_decomposition.h"
#include "qr_eigvals.h"
int main() {
    std::complex<float> a[16] {1, 2, 0, 0.5, 2, -1, 7, 4, 0, 7, 3, 0, 0.5, 4, 0, -2};


    Matrix<std::complex<float>> A(a, 4, 4);
    Matrix<std::complex<float>> Q, R;

    std::pair<Matrix<std::complex<float>>, Matrix<std::complex<float>>> p;
    p = QREigenvalues(A);

    std::cout<<"M1"<<"\n";
    Matrix<std::complex<float>> M1 = p.first;
    for (size_t i = 0; i < M1.VertDim(); i++) {
        for (size_t j = 0; j < M1.HorizDim(); j++) {
            std::cout << M1(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';

    Matrix<std::complex<float>> M2 = p.second;
    for (size_t i = 0; i < M2.VertDim(); i++) {
        for (size_t j = 0; j < M2.HorizDim(); j++) {
            std::cout << M2(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';

}
