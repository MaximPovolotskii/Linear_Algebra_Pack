//
// Created by Максим on 17.11.2021.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <vector>
#include <exception>
#include "matrix_multiplication.h"

const double EPS = 1E-9;


class BadMatrixDimension: public std::exception {
public:
    const char* what() const throw() override {
        return "Dimensions do not converge";
    }
};

class NotSquareMatrix: public BadMatrixDimension {
public:
    const char* what() const throw() override {
        return "Matrix is not square";
    }
};

class IndexOutOfMatrix: public std::exception {
    const char* what() const throw() override {
        return "Index out of range";
    }
};

template<class T>
const T& min(const T& a, const T& b)
{
    return (b < a) ? b : a;
}

enum Matrix_types
{
    NULLMATRIX,  ///нулевая
    IDENTITY, ///единичная
    NONE  ///плевать
};

template <class T>
class Matrix {
private:
    T* coord;
    size_t vert_dim = 0;
    size_t horiz_dim = 0;
public:
    Matrix() = default;
    Matrix(size_t vert_dim_, size_t horiz_dim_, int matrix_type = NONE): coord(new T[vert_dim_*horiz_dim_]) {
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
        if (matrix_type == NULLMATRIX) {
            for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
                coord[i] = 0;
            }
        }
        if (matrix_type == IDENTITY) {
            for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
                coord[i] = 0;
            }
            for (size_t i = 0; i < min(vert_dim, horiz_dim); i++) {
                coord[i*horiz_dim + i] = 1;
            }
        }
    }


    explicit Matrix(T* v, size_t vert_dim_, size_t horiz_dim_): coord(new T[vert_dim_ * horiz_dim_]) {
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
        for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
            coord[i] = v[i];
        }
    }

    Matrix(Matrix<T> &&rhs) noexcept  {
        coord = rhs.coord;
        rhs.coord = nullptr;
        std::swap(vert_dim, rhs.vert_dim);
        std::swap(horiz_dim, rhs.horiz_dim);
    }

    Matrix(const Matrix<T> &rhs) {
        vert_dim = rhs.vert_dim;
        horiz_dim = rhs.horiz_dim;
        coord = new T[rhs.vert_dim * rhs.horiz_dim];
        for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
            coord[i] = rhs.coord[i];
        }
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

    ~Matrix(){
        if (coord != nullptr)   {delete[] coord;}
    }

    std::vector<size_t> Shape() const {
        std::vector<size_t> shape({vert_dim, horiz_dim});
        return shape;
    }

    size_t HorizDim() const{
        return horiz_dim;
    }

    size_t VertDim() const{
        return vert_dim;
    }

    Matrix operator+(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();
        }
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++){
            res[i] = coord[i] + rm.coord[i];
        }
        return Matrix<T>(std::move(res), vert_dim, horiz_dim);
    }

    Matrix operator-(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();
        }
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++){
            res[i] = coord[i] - rm.coord[i];
        }
        return Matrix<T>(std::move(res), vert_dim, horiz_dim);
    }

    Matrix operator*(const T& val) const {
        T* res;
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++) {
            res[i] = coord[i] * val;
        }
        return Matrix<T>(std::move(res), vert_dim, horiz_dim);
    }

    bool operator == (const Matrix<T>& rm) const {
        if (vert_dim != rm.vert_dim || horiz_dim != rm.horiz_dim) {
            return false;
        }
        size_t mat_size = vert_dim * horiz_dim;
        for (size_t i = 0; i < mat_size; i++){
            if (coord[i] != rm.coord[i]) {return false;}
        }
        return true;
    }

    bool operator != (const Matrix<T>& rm) const {
        return !(this == rm);
    }

    T& operator()(const size_t i, const size_t j) const{
        if (i < vert_dim && j < horiz_dim) {
            return coord[i * horiz_dim + j];
        }
        else {
            throw IndexOutOfMatrix();
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

    T Determinant() const {
        if (VertDim() != HorizDim()) {
            throw NotSquareMatrix();
        }
        T det_calc = T(1);
        std::vector<size_t> v;
        Matrix<T> B = StraightRun(*this, Matrix<T>(horiz_dim, 1), det_calc, v).first;
        return det_calc;
    }

    Matrix<T> Inverse() const {
        if (VertDim() != HorizDim()) {
            throw NotSquareMatrix();
        }
        T det_calc = T(0);
        std::pair<Matrix<T>, Matrix<T>> p;
        std::vector<size_t> v;
        p = StraightRun(*this, Matrix<T>(horiz_dim, horiz_dim, IDENTITY), det_calc, v);
        p = ReverseRun(p.first, p.second);
        return p.second;
    };


    size_t Rank() const {
        size_t rank = 0;
        T det_calc = T(0);
        std::vector<size_t> v;
        Matrix<T> B = StraightRun(*this, Matrix<T>(vert_dim, 1), det_calc, v).first;
        for (size_t i = 0; i < min(horiz_dim, vert_dim); i++) {
            if (abs(B(i, i)) > EPS) {
                rank++;
            }
        }
        return rank;
    }

};

#endif //LINALG_MATRIX_H