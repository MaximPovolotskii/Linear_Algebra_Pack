//
// Created by Максим on 17.11.2021.
//

#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <vector>
#include <exception>

class BadMatrixDimension: public std::exception {
public:
    const char* what() const throw() {
        return "Dimensions do not converge";
    }
};

template <class T>
class Matrix {
private:
    std::vector<T> coord{};
    size_t vert_dim = 0;
    size_t horiz_dim = 0;
public:
    Matrix() = default;

    explicit Matrix(const std::vector<T> &v, const size_t &vert_dim_, const size_t &horiz_dim_) {
        coord = v;
        vert_dim = vert_dim_;
        horiz_dim = horiz_dim_;
    }

    explicit Matrix(const std::vector<std::vector<T>> mat) {
        vert_dim = mat.size();
        //Here we need a check whether all horizontal dimensions are correct or not
        //(they must be equal and non-zero)
        //IDK whether to throw an exception or not
        for (size_t i = 0; i < vert_dim; i++) {
            for (auto it = mat.begin(); it != mat.end(); it++) {
                coord.push_back(*it);
            }
        }
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


    std::vector<T> Shape() const {
        std::vector<size_t> shape({vert_dim, horiz_dim});
        return shape;
    }

    Matrix operator+(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();     /// Think about exceptions
        }
        std::vector<T> res;
        for (size_t i = 0; i < coord.size(); i++){
            res.push_back(coord[i] + rm.coord[i]);
        }
        return Matrix<T>(std::move(res));
    }

    Matrix operator-(const Matrix<T>& rm) const {
        if (Shape() != rm.Shape())  {
            throw BadMatrixDimension();     /// Think about exceptions
        }
        std::vector<T> res;
        for (size_t i = 0; i < coord.size(); i++){
            res.push_back(coord[i] - rm.coord[i]);
        }
        return Matrix<T>(std::move(res));
    }

    Matrix operator*(const T& val) const {
        std::vector<T> res;
        for (size_t i = 0; i < coord.size(); i++) {
            res.push_back(coord[i] * val);
        }
        return Matrix<T>(std::move(res));
    }

    bool operator== (const Matrix<T>& rm) const {
        return (coord == rm.coord && vert_dim == rm.vert_dim && horiz_dim == rm.horiz_dim);
    }
    bool operator!= (const Matrix<T>& rm) const {
        return (coord != rm.coord || vert_dim != rm.vert_dim || horiz_dim != rm.horiz_dim);
    }

    T* operator[](size_t n) const { // So that mat[n][m] gives m-1-th element of mat[n]
        if (n < vert_dim) {
            return coord.begin() + (n * horiz_dim);
        }
    }

    Matrix operator* (const Matrix<T>& rm) const {
        if (horiz_dim != rm.vert_dim) {
            throw BadMatrixDimension();
        }

        std::vector<T> res;
        T buffer; //T must be 0-compatible time

        for (size_t i = 0; i < vert_dim; i++) {
            for (size_t j = 0; j < rm.horiz_dim; j++) {
                buffer = 0;
                for (size_t k = 0; k < horiz_dim; k++) {
                    buffer += coord[horiz_dim * i + k] * rm[k][j];
                }
                res.push_back(buffer);
            }
        }

        return Matrix(res, vert_dim, rm.horiz_dim);
    }

};



#endif //LINALG_MATRIX_H
