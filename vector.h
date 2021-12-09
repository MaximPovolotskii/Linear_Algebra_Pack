
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
    Vector();
    Vector(const Matrix<T> & mat);
    Vector(Matrix<T> && mat);
    Vector(const Vector<T> & vect): Matrix<T>(vect) {}
    Vector(Vector<T> && vect): Matrix<T>(vect) {}
    Vector &operator=(Matrix<T> &&rhs) noexcept;
    Vector &operator=(Vector<T> &&rhs) noexcept;
    Vector &operator=(const Matrix<T> &rhs);
    Vector &operator=(const Vector<T> &rhs);
    explicit Vector(size_t n, int type = NONE);
    Vector(T* v, size_t n);
    Vector(std::vector<T> v, size_t n);
    const T& operator() (size_t i) const;
    T& operator() (size_t i);
};

template<class T>
Vector<T>::Vector() : Matrix<T>() {}

template<class T>
Vector<T>::Vector(const Matrix<T> &mat): Matrix<T>(mat){
    if (mat.HorizDim() != 1) {
        std::cerr<<"In Vector(const Matrix<T> & ) ";
        throw BadMatrixDimension();}
}

template<class T>
Vector<T>::Vector(Matrix<T> &&mat): Matrix<T>(mat){
    if (mat.HorizDim() != 1) {
        std::cerr<<"In Vector(Matrix<T> && ) ";
        throw BadMatrixDimension();}
}

template<class T>
Vector<T> &Vector<T>::operator=(Matrix<T> &&rhs) noexcept {
    if (rhs.HorizDim() != 1) {
        std::cerr<<"In Vector &operator=(Matrix<T> &&)   ";
        throw BadMatrixDimension();
    }
    if (this == &rhs) return *this;
    Vector<T> mat(std::move(rhs));
    std::swap(this->Matrix<T>::transposed, mat.transposed);
    std::swap(this->Matrix<T>::coord, mat.coord);
    std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
    std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Vector<T> &Vector<T>::operator=(Vector<T> &&rhs) noexcept {
    if (this == &rhs) return *this;
    Vector<T> mat(std::move(rhs));
    std::swap(this->Matrix<T>::transposed, mat.transposed);
    std::swap(this->Matrix<T>::coord, mat.coord);
    std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
    std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Vector<T> &Vector<T>::operator=(const Matrix<T> &rhs) {
    if (rhs.HorizDim() != 1) {
        std::cerr << "In Vector &operator=(const Matrix<T> &)   ";
        throw BadMatrixDimension();
    }
    if (this == &rhs) return *this;
    Vector<T> mat(rhs);
    std::swap(this->Matrix<T>::transposed, mat.transposed);
    std::swap(this->Matrix<T>::coord, mat.coord);
    std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
    std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs) {
    if (this == &rhs) return *this;
    Vector<T> mat(rhs);
    std::swap(this->Matrix<T>::transposed, mat.transposed);
    std::swap(this->Matrix<T>::coord, mat.coord);
    std::swap(this->Matrix<T>::vert_dim, mat.vert_dim);
    std::swap(this->Matrix<T>::horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Vector<T>::Vector(size_t n, int type) : Matrix<T>(n, 1, type) {}

template<class T>
Vector<T>::Vector(T *v, size_t n) : Matrix<T>(v, n, 1) {}

template<class T>
Vector<T>::Vector(std::vector<T> v, size_t n) : Matrix<T>(v, n, 1) {}

template<class T>
const T &Vector<T>::operator()(size_t i) const {
    if (i >= this->VertDim()) {
        std::cerr<<"In Vector::operator()   ";
        throw BadVectorDimension();
    }
    return this->Matrix<T>::operator()(i, 0);
}

template<class T>
T &Vector<T>::operator()(size_t i) {
    if (i >= this->VertDim()) {
        std::cerr<<"In Vector::operator()   ";
        throw BadVectorDimension();
    }
    return this->Matrix<T>::operator()(i, 0);
}

template<class T>
Vector<T> CrossProduct(Vector<T> lv, Vector<T> rv)  {
    if ((lv.VertDim() != 3)||(rv.VertDim()  != 3))  {
        std::cerr<<"In CrossProduct   ";
        throw BadVectorDimension();
    }
    T a0 = lv(1) * rv(2) - lv(2) * rv(1);
    T a1 = lv(2) * rv(0) - lv(0) * rv(2);
    T a2 = lv(0) * rv(1) - lv(1) * rv(0);
    return Vector<T>({a0, a1, a2}, 3);
}


template <class T>
long double Norm(const Vector<T>& v) {
    long double a = 0;
    for (size_t i = 0; i < v.VertDim(); i++){
        a += std::abs(v(i))*std::abs(v(i));
    }
    a = std::sqrt(a);
    return a;
}

template <class T>
Vector<T> Normalize(const Vector<T>& v){
    T b[v.VertDim()];
    long double n = Norm(v);
    if (std::abs(n) > EPS) {
        for (size_t i = 0; i < v.VertDim(); i++) {
            b[i] = v(i) / n;
        }
    }
    else {
        for (size_t i = 0; i < v.VertDim(); i++) {
            b[i] = 0;
        }
    }
    return Vector<T>(b, v.VertDim());
}

#endif //LINALG_VECTOR_H
