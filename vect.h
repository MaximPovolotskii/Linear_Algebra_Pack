#ifndef LINALG_VECT_H
#define LINALG_VECT_H

#include <vector>
#include <exception>

class BadVectorDimension: public std::exception {
public:
    const char* what() const throw() {
        return "Dimensions do not converge";
    }
};

template <class T>
class Vector {
private:
    std::vector<T> coord {};
public:
    Vector()=default;
    explicit Vector(const std::vector<T> &v): coord(v) {}
    Vector(Vector<T>&& rhs) noexcept : coord(std::move(rhs.coord)) {}
    Vector(const Vector<T>& rhs) : coord(rhs.coord) {}
    Vector& operator= (Vector<T>&& rhs)  noexcept;
    Vector& operator= (const Vector<T>& rhs);
    ~Vector() = default;

    size_t Shape() const {
        return coord.size();
    }

    const T operator[](size_t n) const { ///считать с 1 или с 0?
        return coord[n];
    }

    bool operator== (const Vector<T>& rv) const{
        return (coord == rv.coord);
    }
    bool operator!= (const Vector<T>& rv) const{
        return (coord != rv.coord);
    }

    Vector operator*(const T& val) const {
        std::vector<T> res;
        for (size_t i = 0; i < Shape(); i++){
            res.push_back(coord[i]*val);
        }
        return Vector<T>(res);
    }

    Vector operator+(const Vector<T>& rv) const {
        if (Shape() != rv.Shape())  {
            throw BadVectorDimension();     /// Think about exceptions
        }
        std::vector<T> res;
        for (size_t i = 0; i < Shape(); i++){
            res.push_back(coord[i] + rv.coord[i]);
        }
        return Vector<T>(res);
    }

    Vector operator-(const Vector<T>& rv) const {
        if (Shape() != rv.Shape())  {
            throw BadVectorDimension();
        }
        std::vector<T> res;
        for (size_t i = 0; i < Shape(); i++){
            res.push_back(coord[i] - rv.coord[i]);
        }
        return Vector<T>(res);
    }

    T operator*(const Vector<T>& rv) const {
        if (Shape() != rv.Shape())  {
            throw BadVectorDimension();
        }
        T res = T(0);                       ///Can I do that? Discuss
        for (size_t i = 0; i < Shape(); i++){
            res += coord[i] * rv.coord[i];
        }
        return res;
    }

};

template<class T>
Vector<T> CrossProduct(const Vector<T>& lv, const Vector<T>& rv) {
    if ((lv.Shape() != 3)||(rv.Shape() != 3))  {
        throw BadVectorDimension();
    }
    T a0 = lv.coord[1] * rv.coord[2] - lv.coord[2] * rv.coord[1];
    T a1 = lv.coord[2] * rv.coord[0] - lv.coord[0] * rv.coord[2];
    T a2 = lv.coord[0] * rv.coord[1] - lv.coord[1] * rv.coord[2];
    Vector<T> ans({a0, a1, a2});
    return ans;
}


template<class T>
Vector<T> &Vector<T>::operator=(Vector<T> &&rhs) noexcept {
    if (this == &rhs) return *this;
    Vector<T> t(std::move(rhs));
    std::swap(coord, t.coord);
    return *this;
}

template<class T>
Vector<T> &Vector<T>::operator=(const Vector<T> &rhs) {
    if (this == &rhs) return *this;
    Vector<T> t(rhs);
    std::swap(coord, t.coord);
    return *this;
}

#endif //LINALG_VECT_H
