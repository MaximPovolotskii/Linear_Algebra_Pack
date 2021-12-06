
#ifndef LINALG_LINEAR_SYSTEM_SOLVE_H
#define LINALG_LINEAR_SYSTEM_SOLVE_H

#include "vector.h"
#include "matrix.h"
#include "gauss_solve.h"
#include <tuple>
#include <cmath>

template<typename T>
std::tuple<bool, Vector<T>, Matrix<T>> LinearSystemSolve(const Matrix<T> & A, const Vector<T> & b) {
    T det = 0;
    size_t rank = 0;
    std::vector<size_t> v;
    std::pair<Matrix<T>, Matrix<T>> p;
    p = StraightRun(A, b, det, v, rank);
    p = ReverseRun(p.first, p.second, v, rank);
    Matrix<T> B = std::move(p.first);
    Vector<T> remainder = std::move(p.second);
    for (size_t i = rank; i < remainder.VertDim(); i++) {
        if (abs(remainder(i)) > EPS ) {
            return {false, Vector<T>(A.HorizDim(), NULLMATRIX), Matrix<T>(A.HorizDim(), 1, NULLMATRIX)};   ///  NO SOLUTIONS
        }
    }
    Matrix<T> F(A.HorizDim(), A.HorizDim() - rank, NULLMATRIX);
    for (size_t i  = 0; i < rank; i++) {
        for (size_t j = 0; j < F.HorizDim(); j++) {
            F(i, j) = - B(i, j + rank);
        }
    }
    for (size_t i  = 0; i < F.HorizDim(); i++) {
        F(i + rank, i) = 1;
    }
    std::vector<T> new_coord = {};
    for (size_t i = 0; i < F.VertDim(); i++) {
        if (i < remainder.VertDim()) {
            new_coord.push_back(remainder(i));
        }
        else{
            new_coord.push_back(0);
        }
    }
    Vector<T> ans(new_coord, F.VertDim());
    for (size_t i = 0; i < F.VertDim(); i++) {
        if (v[i] != i) {
            size_t j = v[i];
            for (size_t k = 0; k < F.HorizDim(); k++) {
                std::swap(F(i, k), F(j, k));
            }
            std::swap(remainder(i), remainder(j));
            std::swap(v[i], v[j]);
        }
    }


    return {true, ans, F};


}




#endif //LINALG_LINEAR_SYSTEM_SOLVE_H
