#include <iostream>
#include "matrix.h"

int main() {
    int c[4] = {1, 2, 3, 4};
    Matrix<int> a1(c, 2, 2);
    Matrix<int> b1 = a1.SmartMult(a1);
    std::cout << b1(0, 0) << " " << b1(0, 1) << " " << b1(1, 0) << " " << b1(1, 1) << std::endl;

    float a[12] {0, 0, 0, 0, 2, 0, 3, 4, 7, 0, 8, 1};
    float b[3] {0, 1, 7};
    std::vector<size_t> v;
    Matrix<float> A(a, 3, 4);
    float det = 0;
    std::pair<Matrix<float>, Matrix<float>> p;
    size_t rank = 0;
    p = StraightRun(A, Matrix<float>(b, A.VertDim(), 1), det, v, rank);
    Matrix<float> B = p.first;
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
    Matrix<float> C = p.first;
    for (size_t i = 0; i < C.VertDim(); i++) {
        for (size_t j = 0; j < C.HorizDim(); j++) {
            std::cout << C(i, j) << " ";
        }
        std::cout << '\n';
    }

    return 0;
}
