//
// Created by Максим on 17.11.2021.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <vector>
#include <exception>
#include <iostream>
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
    bool transposed = false;
public:
    Matrix() = default;
    Matrix(size_t vert_dim_, size_t horiz_dim_, int matrix_type = NONE): coord(new T[vert_dim_*horiz_dim_]) {

        transposed = false;
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
        transposed = false;
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
        for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
            coord[i] = v[i];
        }
    }

    Matrix (std::vector<T> v, size_t vert_dim_, size_t horiz_dim_):  coord(new T[vert_dim_ * horiz_dim_]){
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
        if (v.size() != vert_dim_*horiz_dim_) {
            throw BadMatrixDimension();
        }
        for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
            coord[i] = v[i];
        }
    }


    Matrix(Matrix<T> &&rhs) noexcept  {
        coord = rhs.coord;
        rhs.coord = nullptr;
        std::swap(transposed, rhs.transposed);
        std::swap(vert_dim, rhs.vert_dim);
        std::swap(horiz_dim, rhs.horiz_dim);
    }

    Matrix(const Matrix<T> &rhs) {
        transposed = rhs.transposed;
        vert_dim = rhs.vert_dim;
        horiz_dim = rhs.horiz_dim;
        coord = new T[rhs.vert_dim * rhs.horiz_dim];
        for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
            coord[i] = rhs.coord[i];
        }
    }

    Matrix &operator=(Matrix<T> &&rhs) noexcept {
        if (this == &rhs) return (*this);
        Matrix<T> mat(std::move(rhs));
        std::swap(transposed, mat.transposed);
        std::swap(coord, mat.coord);
        std::swap(vert_dim, mat.vert_dim);
        std::swap(horiz_dim, mat.horiz_dim);
        return (*this);
    }

    Matrix &operator=(const Matrix<T> &rhs) {
        if (this == &rhs) return (*this);
        Matrix<T> mat(rhs);
        std::swap(transposed, mat.transposed);
        std::swap(coord, mat.coord);
        std::swap(vert_dim, mat.vert_dim);
        std::swap(horiz_dim, mat.horiz_dim);
        return (*this);
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
        for (size_t i = 0; i < vert_dim; ++i) {
            for (size_t j = 0; j < horiz_dim; ++j) {
                res[i*horiz_dim + j] = (*this)(i, j) + rm(i, j);
            }
        }
        return Matrix<T>(std::move(res), vert_dim, horiz_dim);
    }

    Matrix operator-(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();
        }
        T* res;
        for (size_t i = 0; i < vert_dim; ++i) {
            for (size_t j = 0; j < horiz_dim; ++j) {
                res[i*horiz_dim + j] = (*this)(i, j) - rm(i, j);
            }
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
        for (size_t i = 0; i < vert_dim; ++i) {
            for (size_t j = 0; j < horiz_dim; j++) {
                if ((*this)(i, j) != rm(i, j)) { return false; }
            }
        }
        return true;
    }

    bool operator != (const Matrix<T>& rm) const {
        return this != rm;
    }

    T& operator()(const size_t i, const size_t j) const{
        if (i < vert_dim && j < horiz_dim) {
            if (transposed)
                return coord[j * vert_dim + i];
            else
                return coord[i * horiz_dim + j];
        }
        /*else {
            throw IndexOutOfMatrix();
        }*/
    }

    Matrix operator* (Matrix<T>& rm) {
        /*if (horiz_dim != rm.vert_dim) {
            throw BadMatrixDimension();
        }*/

        T* res = new T[rm.horiz_dim * vert_dim];
        T buffer; //T must be 0-compatible type

        for (size_t i = 0; i < vert_dim; i++) {
            for (size_t j = 0; j < rm.horiz_dim; j++) {
                buffer = 0;
                for (size_t k = 0; k < horiz_dim; k++) {
                    buffer += (((*this)))(i, k) * rm(k, j);
                }
                res[i * rm.horiz_dim + j] = buffer;
            }
        }
        return Matrix(res, vert_dim, rm.horiz_dim);
    }

    Matrix<T> Pow(Matrix<T>& rm, int n){///without ln working time
        if (n==0){
            T* res = new T[rm.horiz_dim * rm.vert_dim];
            for (int i=0;i<rm.horiz_dim;i++){
                for (size_t j = 0; j < rm.horiz_dim; j++){
                    if (i==j){res[i * horiz_dim + j] = 1;
                    }
                    else {res[i * horiz_dim + j] = 0;}
                }

            }
            return Matrix(res, rm.vert_dim, rm.horiz_dim);
        }
        else{
        if (vert_dim==horiz_dim){
            T* res = new T[rm.horiz_dim * rm.vert_dim];
            for (size_t i = 0; i < rm.vert_dim; i++) {
                for (size_t j = 0; j < rm.horiz_dim; j++){
                    res[i * horiz_dim + j] = rm(i,j);
                }
            }
            T buffer; //T must be 0-compatible type
            for (size_t l = 1; l-1 <= n-2; l++) {
                for (size_t i = 0; i < vert_dim; i++) {
                    for (size_t j = 0; j < rm.horiz_dim; j++) {
                        buffer = 0;
                        for (size_t k = 0; k < horiz_dim; k++) {
                            buffer += res[i * rm.horiz_dim + k] * rm(k, j);
                        }
                        res[i * rm.horiz_dim + j] = buffer;
                    }
                }

            }
            //if (l!=n-1){rm = Matrix(res, vert_dim, rm.horiz_dim);}
            return Matrix(res, rm.vert_dim, rm.horiz_dim);
            }}
            //return (*this);}
        }
        ///else - here need to be an exeption

    Matrix<T> PowersTwo(Matrix<T> rm, int a){///4 is equal to [0,0,1]
        Matrix<T> a_2[a];
        int b_2[a];
        int i = 0;
        while (a > 0) {
            //PrintMatrixInt(rm.Pow(rm,i));
            a_2[i] = rm.Pow(rm,i)*(a % 2);
            a = a / 2;
            i++;
        }
        Matrix<T> Result = rm;
        for (int k=0;k<i;k++){

            Result = Result * a_2[k];
        };
        return Result;
    }

    Matrix SmartMult(Matrix<T>& rm) { //to fix non-8-divisible cases and transposed cases
        if (horiz_dim != rm.vert_dim) {
            throw BadMatrixDimension();
        }

        size_t M = (vert_dim/8 + 1) * 8;
        size_t K = (horiz_dim/8 + 1) * 8;
        size_t N = (rm.horiz_dim/8 + 1) * 8;

        T* left_mat = ((*this)).Expand_And_HardTranspose(M, K, false);
        T* right_mat = rm.Expand_And_HardTranspose(K, N, true);

        T* res = new T[M * N];
        T* res_norm = new T[vert_dim * rm.horiz_dim];
        gemm_v2(M, K, N, left_mat, right_mat, res);

        for (size_t i = 0; i < vert_dim; ++i) {
            size_t cur_raw = i * rm.horiz_dim;
            size_t cur_8raw = i * N;
            for (size_t j = 0; j < rm.horiz_dim; ++j) {
                res_norm[cur_raw + j] = res[cur_8raw + j];
            }
        }
        return Matrix(res_norm, vert_dim, rm.horiz_dim);
    }

    T Determinant() const {
        if (VertDim() != HorizDim()) {
            throw NotSquareMatrix();
        }
        T det_calc = T(1);
        std::vector<size_t> v;
        size_t rank = 0;
        Matrix<T> B = StraightRun((*this), Matrix<T>(horiz_dim, 1), det_calc, v, rank).first;
        return det_calc;
    }

    Matrix<T> Inverse() const {
        if (VertDim() != HorizDim()) {
            throw NotSquareMatrix();
        }
        T det_calc = T(0);
        size_t rank = 0;
        std::pair<Matrix<T>, Matrix<T>> p;
        std::vector<size_t> v;

        p = StraightRun((*this), Matrix<T>(horiz_dim, horiz_dim, IDENTITY), det_calc, v, rank);
        p = ReverseRun(p.first, p.second, v, rank);
        return p.second;
    };


    size_t Rank() const {
        size_t rank = 0;
        T det_calc = T(0);
        std::vector<size_t> v;
        Matrix<T> B = StraightRun((*this), Matrix<T>(vert_dim, 1), det_calc, v, rank).first;
        return rank;
    }

    void Transpose() {
        transposed = !transposed;
        std::swap(vert_dim, horiz_dim);
    }

    T* ExpandAndHardTranspose(size_t M, size_t K, bool flip) { ///here need to be CamelCase
        T* left_mat;
        left_mat = new T[M * K];
        if (!flip) {
            for (size_t i = 0; i < M; ++i) {
                size_t cur_raw = i * K;
                if (i < vert_dim) {
                    for (size_t j = 0; j < K; ++j) {
                        if (j < horiz_dim)
                            left_mat[cur_raw + j] = ((*this))(i, j);
                        else
                            left_mat[cur_raw + j] = 0;
                    }
                } else {
                    for (size_t j = 0; j < K; ++j) {
                        left_mat[cur_raw + j] = 0;
                    }
                }
            }
        }
        else {
            for (size_t j = 0; j < K; ++j) {
                size_t cur_col = j * M;
                if (j < horiz_dim) {
                    for (size_t i = 0; i < M; ++i) {
                        if (i < vert_dim)
                            left_mat[cur_col + i] = ((*this))(i, j);
                        else
                            left_mat[cur_col + i] = 0;
                    }
                } else {
                    for (size_t i = 0; i < M; ++i) {
                        left_mat[cur_col + i] = 0;
                    }
                }
            }
        }
        return left_mat;
    }
};


#endif //LINALG_MATRIX_H