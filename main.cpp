#include <iostream>
#include "matrix.h"
#include "vector.h"
#include <complex>
#include <cmath>
#include "qr_decomposition.h"
#include "gauss_run.h"
int main() {
    float a[35] {0, 0, 0, 0, 0, 1, 2, 4, 8, -3, 4, 5, -3, 0, 0, 6, 3, 1, -1, 0, 4, -7, 1, 2, -5, 2, 0, 1, 0, 0, 4, 0, 2, 0, 0};
    float b[7] {0, 12, 7, 5, -1, 1, 2};

    std::vector<size_t> v;
    Matrix<float> A(a, 7, 5);
    Vector<float> s(b, 7);

    float det = 0;
    size_t rank = 0;


    std::pair<Matrix<float>, Matrix<float>> p;
    p = StraightRun(A, s, det, v, rank);
    Matrix<float> B = (p.first);
    Vector<float> remainder = (p.second);
    for (size_t i = 0; i < B.VertDim(); i++) {
        for (size_t j = 0; j < B.HorizDim(); j++) {
            std::cout << B(i, j) << " ";
        }
        std::cout << '\n';
    }
    std::cout << '\n';
    std::cout << '\n';

    p = ReverseRun(p.first, p.second, v, rank);
    B = std::move(p.first);
    remainder = std::move(p.second);
    for (size_t i = 0; i < B.VertDim(); i++) {
        for (size_t j = 0; j < B.HorizDim(); j++) {
            std::cout << B(i, j) << " ";
        }
        std::cout << '\n';
    }



    return 0;

}
