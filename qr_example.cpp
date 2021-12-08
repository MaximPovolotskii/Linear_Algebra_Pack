#include <iostream>
#include <vector>
#include "vector.h"
#include "matrix.h"
#include "gauss_run.h"
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
    std::cout<<"Eigenvectors test"<<"\n";
    for (auto a : v) {
        Matrix<long double> A(a, 4, 4);

        Matrix<long double> Q(a, 4, 4);
        Matrix<long double> R(a, 4, 4);

        auto pqr = (QREigenvalues(A));
        std::vector<long double> evr = pqr.first;
        std::vector<std::complex<long double>> evi = pqr.second;

        std::cout << "real eigenvalues"<<"\n";
        for (long double i : evr) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        std::cout << "complex eigenvalues"<<"\n";
        for (std::complex<long double> i : evi) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        for (long double e : evr) {
            auto p = FindEigenvectors(A, e);
            std::cout<<"eigenvalue "<<e<<"\n"<<"number of eigenvectors "<<p.first<<"\n";
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
            std::cout<<"eigenvalue"<<e<<"\n"<<"number of eigenvectors  "<<p.first<<"\n";
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
    std::cout << '\n' << "\n";

    std::cout<<"Generalized eigenvectors test"<<"\n";
    std::vector<long double*> u;
    long double c0[9] {0, -3, 3, -1, -4, 6, 0, -2, 2};
    u.push_back(c0);
    for (auto c : u) {
        Matrix<long double> C(c, 3, 3);

        auto pqr = (QREigenvalues(C));
        std::vector<long double> evr = pqr.first;
        std::vector<std::complex<long double>> evi = pqr.second;

        std::cout << "real eigenvalues"<<"\n";
        for (long double i : evr) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        std::cout << "complex eigenvalues"<<"\n";
        for (std::complex<long double> i : evi) {
            std::cout << i <<" ";
        }
        std::cout << '\n';

        long double e = evr[0];
        auto p0 = FindEigenvectors(C, e);
        std::cout<<"eigenvalue "<<e<<"\n"<<"number of eigenvectors "<<p0.first<<"\n";
        auto B = p0.second;
        for (size_t i = 0; i < B.VertDim(); i++) {
            for (size_t j = 0; j < B.HorizDim(); j++) {
                std::cout << B(i, j) << " ";
            }
            std::cout << '\n';
        }
        std::cout << '\n';

        Vector<long double> V0(B);

        auto t1 = FindGeneralizedEigenvectors(C, e, V0);
        std::cout<<"eigenvalue "<<e<<"\n"<<"number of joined eigenvectors "<<std::get<0>(t1)<<"\n";
        auto V1 = std::get<1>(t1);
        for (size_t i = 0; i < V1.VertDim(); i++) {
            std::cout << V1(i) << " ";
            std::cout << '\n';
        }
        std::cout << '\n';


    }
    return 0;
}