#include <iostream>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "gauss_solve.h"
#include "qr_decomposition.h"
int main() {
    float a[15] {12, -51, 4, 6, 167, -68, -4, 24, -41, -1, 1, 0, 2, 0, 3};
    float b[15] {3, 2, -1, 5, -4, 24, -41, -1, 1, 0, 2, 0, 3};
    Vector<float> V1(a, 3);
    Vector<float> V2(b, 3);


    Matrix<float> A(a, 5, 3);
    Matrix<float> Q, R;
    Householder(A, R, Q);
    std::cout<<"R"<<"\n";
    ShowMatrix(R);
    std::cout<<"Q"<<"\n";
    ShowMatrix(Q);
}
