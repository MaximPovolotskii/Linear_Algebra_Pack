#include <iostream>
#include "matrix.h"
#include "gauss_solve.h"
#include <complex>
#include <cmath>

int main() {
    std::complex<float> a[12] {0, 0, 0, 0, 2, 0, 3, 4, 7, 0, 8, 1};
    std::complex<float> b[3] {0, 1, 7};
    std::vector<size_t> v;
    Matrix<std::complex<float>> A(a, 3, 4);
    std::complex<float> det = 0;
    std::pair<Matrix<std::complex<float>>, Matrix<std::complex<float>>> p;
    size_t rank = 0;
    p = StraightRun(A, Matrix<std::complex<float>>(b, A.VertDim(), 1), det, v, rank);
    Matrix<std::complex<float>> B = p.first;
    for (size_t i = 0; i < B.VertDim(); i++) {
        for (size_t j = 0; j < B.HorizDim(); j++) {
            std::cout << B(i, j) << " ";
        }
        std::cout << '\n';
    }
    try {
        p = ReverseRun(p.first, p.second, v, rank);
    }
    catch(std::exception &e) {
        std::cerr<<e.what();
    }
    std::cout<< A.Rank();

    std::cout << '\n';
    Matrix<std::complex<float>> C = p.first;
    for (size_t i = 0; i < C.VertDim(); i++) {
        for (size_t j = 0; j < C.HorizDim(); j++) {
            std::cout << C(i, j) << " ";
        }
        std::cout << '\n';
    }

    std::cout << '\n';

    C.Transpose();
    for (size_t i = 0; i < C.VertDim(); i++) {
        for (size_t j = 0; j < C.HorizDim(); j++) {
            std::cout << C(i, j) << " ";
        }
        std::cout << '\n';
    }

    std::cout << '\n';

    std::complex<float> c[4] = {1, 2, 3, 4};
    Matrix<std::complex<float>> mat_a(c, 2, 2);
    Matrix<std::complex<float>> mat_b = mat_a.SmartMult(mat_a);
    std::cout << mat_b(0, 0) << " " << mat_b(0, 1) << "\n" << mat_b(1, 0) << " " << mat_b(1, 1) << std::endl;
    return 0;

    return 0;
}
