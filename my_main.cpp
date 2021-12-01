#include <iostream>
#include "gauss_solve.h"
#include "matrix.h"
#include "vector.h"
#include "print.h"


int main() {
    int n=3, m=4, N = 12;
    int max;
    int a[12] {1000, 2530, 0, 0, 2, 260, 3, 4, 70, 0, 8, 1};
    /*for (int i = 1; i < N; i++) {
        if (a[i] > max) {
            max = a[i];
        }
    }*/

    Matrix<int> A(a, n, m);
    max = GiveAbsInt(A);
    PrintMatrixInt(A);
    std::cout<<GiveTrace(A)<<std::endl;

    float b[9] {0,0,1,0,0,0,0,0,1};
    //Matrix<int> B(b, n, n);
    //PrintMatrixInt(B);
    //std::cout<<GiveTrace(B)<<std::endl;

    Vector<float> v(b,n*3);
    PrintVectorFloat(Norm(v));
}
