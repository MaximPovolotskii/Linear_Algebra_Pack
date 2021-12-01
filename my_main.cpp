#include <iostream>
#include "gauss_solve.h"
#include "matrix.h"
#include "vect.h"
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
    max = give_abs_int(A, n, m);
    print_matrix_int(A,n,m);
    std::cout<<give_trace(A,n,m)<<std::endl;

    int b[12] {1000, 2530, 0, 0, 2, 260, 3, 4, 70};
    Matrix<int> B(b, n, n);
    print_matrix_int(B,n,n);
    std::cout<<give_trace(B,n,n)<<std::endl;

}