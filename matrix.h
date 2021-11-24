//
// Created by Максим on 17.11.2021.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <vector>
#include <exception>
#include "matrix_multiplication.h"

class BadMatrixDimension: public std::exception {
public:
    const char* what() const throw() {
        return "Dimensions do not converge";
    }
};

template <class T>
class Matrix {
private:
    T* coord{};
    size_t vert_dim = 0;
    size_t horiz_dim = 0;
public:
    Matrix() = default;

    explicit Matrix(T* v, const size_t &vert_dim_, const size_t &horiz_dim_) {
        coord = v;
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
    }

    Matrix(Matrix<T> &&rhs) noexcept : coord(std::move(rhs.coord)) {
        std::swap(vert_dim, rhs.vert_dim);
        std::swap(horiz_dim, rhs.horiz_dim);
    }

    Matrix(const Matrix<T> &rhs) : coord(rhs.coord) {
        vert_dim = rhs.vert_dim;
        horiz_dim = rhs.horiz_dim;
    }

    Matrix &operator=(Matrix<T> &&rhs) noexcept {
        if (this == &rhs) return *this;
        Matrix<T> mat(std::move(rhs));
        std::swap(coord, mat.coord);
        std::swap(vert_dim, mat.vert_dim);
        std::swap(horiz_dim, mat.horiz_dim);
        return *this;
    }

    Matrix &operator=(const Matrix<T> &rhs) {
        if (this == &rhs) return *this;
        Matrix<T> mat(rhs);
        std::swap(coord, mat.coord);
        std::swap(vert_dim, mat.vert_dim);
        std::swap(horiz_dim, mat.horiz_dim);
        return *this;
    }

    ~Matrix() = default;


    T* Shape() const {
        std::vector<size_t> shape({vert_dim, horiz_dim});
        return shape;
    }

    Matrix operator+(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();     /// Think about exceptions
        }
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++){
            res[i] = coord[i] + rm.coord[i];
        }
        return Matrix<T>(std::move(res));
    }

    Matrix operator-(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();     /// Think about exceptions
        }
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++){
            res[i] = coord[i] - rm.coord[i];
        }
        return Matrix<T>(std::move(res));
    }

    Matrix operator*(const T& val) const {
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++) {
            res[i] = coord[i] * val;
        }
        return Matrix<T>(std::move(res));
    }

    bool operator== (const Matrix<T>& rm) const {
        return (coord == rm.coord && vert_dim == rm.vert_dim && horiz_dim == rm.horiz_dim);
    }
    bool operator!= (const Matrix<T>& rm) const {
        return (coord != rm.coord && vert_dim == rm.vert_dim && horiz_dim == rm.horiz_dim);
    }

    T& operator()(const size_t i, const size_t j) {
        if (i < vert_dim && j < horiz_dim) {
            return coord[i * horiz_dim + j];
        }
    }

    Matrix operator* (Matrix<T>& rm) {
        if (horiz_dim != rm.vert_dim) {
            throw BadMatrixDimension();
        }

        T* res = new T[rm.horiz_dim * vert_dim];
        T buffer; //T must be 0-compatible type

        for (size_t i = 0; i < vert_dim; i++) {
            for (size_t j = 0; j < rm.horiz_dim; j++) {
                buffer = 0;
                for (size_t k = 0; k < horiz_dim; k++) {
                    buffer += coord[horiz_dim * i + k] * rm(k, j);
                }
                res[i * rm.horiz_dim + j] = buffer;
            }
        }
        return Matrix(res, vert_dim, rm.horiz_dim);
    }

    Matrix SmartMult(Matrix<T>& rm) {
        if (horiz_dim != rm.vert_dim) {
            throw BadMatrixDimension();
        }

        T* res = new T[rm.horiz_dim * vert_dim];
        gemm_v2(vert_dim, horiz_dim, rm.horiz_dim, coord, rm.coord, res);

        return Matrix(res, vert_dim, rm.horiz_dim);
    }
};





#endif //LINALG_MATRIX_H