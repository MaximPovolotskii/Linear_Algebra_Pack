#ifndef LINALG_QR_DECOMPOSITION_H
#define LINALG_QR_DECOMPOSITION_H

#include "matrix.h"
#include <vector>

template<typename T>
Matrix<T> ComputeHouseholderFactor(const Vector<T>& v)
{
    size_t n = v.VertDim();
    Matrix<T> res(n, n, NULLMATRIX);
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            res(i,j) = -2 *  v(i) * v(j);   ///conj
    for (size_t i = 0; i < n; i++)
        res(i,i) += 1;
    return res;
}

template<typename T>
void ExtractMinor(const Matrix<T>& mat, size_t d, Matrix<T> & res) {
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

template<typename T>
void Householder(const Matrix<T>& mat, Matrix<T>& R, Matrix<T>& Q) {
    size_t vert_dim = mat.VertDim();
    size_t horiz_dim = mat.HorizDim();
    std::vector<Matrix<T>> qv(vert_dim);

    Matrix<T> z(mat);
    Matrix<T> z1(vert_dim, horiz_dim);
    Vector<T> e(vert_dim), x(vert_dim);
    for (size_t k = 0; k < horiz_dim && k < vert_dim - 1; k++) {
        ExtractMinor(z, k, z1);
        ExtractColumn(z1, x, k);
        T a = Norm(x);
        if (mat(k,k) > 0) {a = -a;}
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
