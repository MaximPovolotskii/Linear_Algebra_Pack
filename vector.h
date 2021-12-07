
#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include "matrix.h"
#include <cmath>
#include <complex>

class BadVectorDimension: public BadMatrixDimension {
public:
    const char* what() const throw() override {
        return "Dimensions of vectors do not converge";
    }
};

template<class T>
class Vector : public Matrix<T> {
public:
    Vector() : Matrix<T>() {};
    Vector(const Matrix<T> & mat): Matrix<T>(mat){
        if (mat.HorizDim() != 1) {throw BadMatrixDimension();}
    }
    Vector(Matrix<T> && mat): Matrix<T>(mat){
        if (mat.HorizDim() != 1) {throw BadMatrixDimension();}
    }
    Vector(const Vector<T> & vect): Matrix<T>(vect) {}
    Vector(Vector<T> && vect): Matrix<T>(vect) {}
    Vector &operator=(Matrix<T> &&rhs) noexcept {
        if (this == &rhs) return *this;
        Vector<T> mat(std::move(rhs));
        std::swap(this->Matrix<T>::transposed, mat.transposed);
        std::swap(this->Matrix<T>::coord, mat.coord);
        std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
        std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
        return *this;
    }
    Vector &operator=(Vector<T> &&rhs) noexcept {
        if (this == &rhs) return *this;
        Vector<T> mat(std::move(rhs));
        std::swap(this->Matrix<T>::transposed, mat.transposed);
        std::swap(this->Matrix<T>::coord, mat.coord);
        std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
        std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
        return *this;
    }
    Vector &operator=(const Matrix<T> &rhs) {
        if (this == &rhs) return *this;
        Vector<T> mat(rhs);
        std::swap(this->Matrix<T>::transposed, mat.transposed);
        std::swap(this->Matrix<T>::coord, mat.coord);
        std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
        std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
        return *this;
    }
    Vector &operator=(const Vector<T> &rhs) {
        if (this == &rhs) return *this;
        Vector<T> mat(rhs);
        std::swap(this->Matrix<T>::transposed, mat.transposed);
        std::swap(this->Matrix<T>::coord, mat.coord);
        std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
        std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
        return *this;
    }
    explicit Vector(size_t n, int type = NONE) : Matrix<T>(n, 1, type) {}
    Vector(T* v, size_t n) : Matrix<T>(v, n, 1) {}
    Vector(std::vector<T> v, size_t n) : Matrix<T>(v, n, 1) {}
    const T& operator() (size_t i) const {
        if (i >= this->VertDim()) {
            std::cerr<<"Vector::operator()";
            throw BadVectorDimension();
        }
        return this->Matrix<T>::operator()(i, 0);
    }

    T& operator() (size_t i) {
        if (i >= this->VertDim()) {
            std::cerr<<"Vector::operator()";
            throw BadVectorDimension();
        }
        return this->Matrix<T>::operator()(i, 0);
    }
};

template<class T>
Vector<T> CrossProduct(Vector<T> lv, Vector<T> rv)  {
    if ((lv.VertDim() != 3)||(rv.VertDim()  != 3))  {
        throw BadVectorDimension();
    }
    T a0 = lv(1) * rv(2) - lv(2) * rv(1);
    T a1 = lv(2) * rv(0) - lv(0) * rv(2);
    T a2 = lv(0) * rv(1) - lv(1) * rv(0);
    return Vector<T>({a0, a1, a2}, 3);
}


template <class T>
T Norm(const Vector<T>& v) {
    T a = 0;
    for (size_t i = 0; i < v.VertDim(); i++){
        a += std::abs(v(i))*std::abs(v(i));
    }
    a = std::sqrt(a);
    return a;
}

template <class T>
Vector<T> Normalize(const Vector<T>& v){
    T b[v.VertDim()];
    T n = Norm(v);
    if (std::abs(n) > EPS) {
        for (size_t i = 0; i < v.VertDim(); i++) {
            b[i] = v(i) / n;
        }
    }
    return Vector<T>(b, v.VertDim());
}

#endif //LINALG_VECTOR_H
