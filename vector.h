
#ifndef LINALG_VECTOR_H
#define LINALG_VECTOR_H

#include "matrix.h"



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
    explicit Vector(size_t n, int type = NONE) : Matrix<T>(n, 1, type) {}
    Vector(T* v, size_t n) : Matrix<T>(v, n, 1) {}
    Vector(std::vector<T> v, size_t n) : Matrix<T>(v, n, 1) {}
    T& operator() (size_t i) {
        /*if (i >= this->VertDim()) {
            throw BadVectorDimension();
        }*/
        return this->Matrix<T>::operator()(i, 0);
    }

};

template<class T>
Vector<T> CrossProduct(Vector<T> lv, Vector<T> rv)  {
    /*if ((lv.VertDim() != 3)||(rv.VertDim()  != 3))  {
        throw BadVectorDimension();
    }*/
    T a0 = lv(1) * rv(2) - lv(2) * rv(1);
    T a1 = lv(2) * rv(0) - lv(0) * rv(2);
    T a2 = lv(0) * rv(1) - lv(1) * rv(0);
    return Vector<T>({a0, a1, a2}, 3);
}


template<class T>
Vector<T> Norm(Vector<T> v){
    float a=0;
    for (uint32_t i =0; i<v.VertDim();i++){
        a += v(i)*v(i);
    }
    a = sqrt(a);
    T b[v.VertDim()];
    for (uint32_t i =0; i<v.VertDim();i++){
        b[i] = v(i)/a;
    }
    return Vector<T>(b,v.VertDim());
}


#endif //LINALG_VECTOR_H
