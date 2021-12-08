#include <iostream>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "gauss_solve.h"
#include "qr_decomposition.h"
#include "qr_eigvals.h"
#include <complex>
#include <tuple>


using namespace std::complex_literals;
int main() {


    std::vector<long double*> v;
    long double b0[16] {1, 2, 0, 1, 0, -2, 7, 4, 4, 0, 3, 0, -6, 0, 6, 8};
    long double b1[16] {6, 9, 0, 4, 4, 9, -7, 0, 0, 0, 0, 1, 3, -12, 0, 0};
    v.push_back(b0);
    v.push_back(b1);
    for (auto a : v) {
        Matrix<long double> A(a, 4, 4);

        Matrix<long double> Q(a, 4, 4);
        Matrix<long double> R(a, 4, 4);


        auto pqr = (QREigenvalues(A));
        std::vector<long double> evr = pqr.first;
        std::vector<std::complex<long double>> evi = pqr.second;

        std::cout << "real eigenvals"<<"\n";
        for (long double i : evr) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        std::cout << "complex eigenvals"<<"\n";
        for (std::complex<long double> i : evi) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        for (long double e : evr) {
            auto p = FindEigenvectors(A, e);
            std::cout<<"eigvalue "<<e<<"\n"<<"number of eigvectors "<<p.first<<"\n";
            auto B = p.second;
            for (size_t i = 0; i < B.VertDim(); i++) {
                for (size_t j = 0; j < B.HorizDim(); j++) {
                    std::cout << B(i, j) << " ";
                }
                std::cout << '\n';
            }
            std::cout << '\n';
        }

        for (std::complex<long double> e : evi) {
            auto p = FindEigenvectors(A, e);
            std::cout<<"eigvalue"<<e<<"\n"<<"number of eigvectors  "<<p.first<<"\n";
            auto B = p.second;
            for (size_t i = 0; i < B.VertDim(); i++) {
                for (size_t j = 0; j < B.HorizDim(); j++) {
                    std::cout << B(i, j) << " ";
                }
                std::cout << '\n';
            }
            std::cout << '\n';
        }
    }


    return 0;
}