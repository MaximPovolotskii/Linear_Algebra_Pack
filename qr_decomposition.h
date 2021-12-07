#ifndef LINALG_QR_DECOMPOSITION_H
#define LINALG_QR_DECOMPOSITION_H

#include "matrix.h"
#include <vector>
#include <complex>
#include <type_traits>

template<typename CT>
Matrix<CT> ComputeHouseholderFactor(const Vector<CT>& v)
{
    size_t n = v.VertDim();
    Matrix<CT> res(n, n, NULLMATRIX);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            res(i,j) = v(i) * std::conj(v(j)) * static_cast<CT>(-2.);   ///conj
    for (size_t i = 0; i < n; i++)
        res(i,i) += 1;
    return res;
}

template<typename CT>
void ExtractMinor(const Matrix<CT>& mat, size_t d, Matrix<CT> & res) {
    for (size_t i = 0; i < mat.VertDim(); i++) {
        for (size_t j = 0; j < mat.HorizDim(); j++)
            res(i,j) = 0;
    }
    for (size_t i = 0; i < d; i++)
            res(i,i) = 1;
    for (size_t i = d; i < mat.VertDim(); i++) {
        for (size_t j = d; j < mat.HorizDim(); j++)
            res(i,j) = mat(i,j);
    }
}

template<typename T>
void ExtractColumn(const Matrix<T> &mat, Vector<T>& v, size_t c){
    for (size_t i = 0; i < mat.VertDim(); i++)
        v(i) = mat(i, c);
}

template<typename CT>
void QRDecomposition(const Matrix<CT>& mat, Matrix<CT>& Q, Matrix<CT>& R) {
    size_t vert_dim = mat.VertDim();
    size_t horiz_dim = mat.HorizDim();
    std::vector<Matrix<CT>> qv(vert_dim);

    Matrix<CT> z(mat);
    Matrix<CT> z1(vert_dim, horiz_dim);
    Vector<CT> e(vert_dim), x(vert_dim);
    for (size_t k = 0; k < horiz_dim && k < vert_dim - 1; k++) {
        ExtractMinor(z, k, z1);
        ExtractColumn(z1, x, k);
        CT a = Norm(x);
        if (std::real(mat(k,k)) > 0) {a = -a;}
        for (size_t i = 0; i < e.VertDim(); i++)
            e(i) = (i == k) ? 1 : 0;
        e = x + e*a;
        e = Normalize(e);
        // qv[k] = I - 2 *e*e^T
        qv[k] = ComputeHouseholderFactor(e);
        z = (qv[k] * z1);

    }
    Q = qv[0];

    for (size_t i = 1; i < horiz_dim && i < vert_dim - 1; i++) {
        Q = qv[i] * Q;
    }
    R = Q*mat;
    Q.Transpose();
}



#endif //LINALG_QR_DECOMPOSITION_H
