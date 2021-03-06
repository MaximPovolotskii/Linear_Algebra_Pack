
#ifndef LINALG_MATRIX_H
#define LINALG_MATRIX_H
#include <vector>
#include <exception>
#include "decl.h"

template <class T>
class Matrix;


#include "matrix_multiplication.h"
#include "gauss_run.h"
#include <complex>
const double EPS = 1E-10;

template<typename T>
struct is_complex : public std::false_type {};

template<typename T>
struct is_complex<std::complex<T>> : public std::true_type {};

enum Matrix_types
{
    NULLMATRIX,  ///нулевая
    IDENTITY, ///единичная
    NONE  ///плевать
};

template <class T>
class Matrix {
protected:
    T* coord = nullptr;
    size_t vert_dim = 0;
    size_t horiz_dim = 0;
    bool transposed = false;
public:
    Matrix() = default;
    Matrix(size_t vert_dim_, size_t horiz_dim_, int matrix_type = NONE);
    explicit Matrix(T* v, size_t vert_dim_, size_t horiz_dim_);
    Matrix (std::vector<T> v, size_t vert_dim_, size_t horiz_dim_);
    Matrix(Matrix<T> &&rhs) noexcept;
    Matrix(const Matrix<T> &rhs);
    virtual Matrix &operator=(Matrix<T> &&rhs) noexcept;
    virtual Matrix &operator=(const Matrix<T> &rhs);
    ~Matrix();
    std::vector<size_t> Shape() const;
    size_t HorizDim() const;
    size_t VertDim() const;
    Matrix operator+(const Matrix<T>& rm) const;
    Matrix operator-(const Matrix<T>& rm) const;
    Matrix operator*(const T& val) const;
    bool operator== (const Matrix<T>& rm) const;
    bool operator!= (const Matrix<T>& rm) const;
    const T& operator()(size_t i, size_t j) const;
    T& operator()(size_t i, size_t j);
    Matrix operator* (const Matrix<T>& rm) const;
    Matrix SmartMult(Matrix<T>& rm);
    T Determinant() const;
    size_t Rank() const;
    void Transpose();
    T* ExpandAndHardTranspose(size_t M, size_t K, bool flip);
    template <typename Q = T>
    std::enable_if_t<is_complex<Q>::value, void> Conjugate();
    template <typename Q = T>
    std::enable_if_t<!is_complex<Q>::value, void> Conjugate();
    Matrix<T> Pow(Matrix<T> rm, int a);///else - here need to be an exception

};

template<class T>
Matrix<T>::Matrix(size_t vert_dim_, size_t horiz_dim_, int matrix_type): coord(new T[vert_dim_*horiz_dim_]) {

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

template<class T>
Matrix<T>::Matrix(T *v, size_t vert_dim_, size_t horiz_dim_): coord(new T[vert_dim_ * horiz_dim_]) {
    transposed = false;
    vert_dim = vert_dim_;
    horiz_dim = horiz_dim_;
    for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
        coord[i] = v[i];
    }
}

template<class T>
Matrix<T>::Matrix(std::vector<T> v, size_t vert_dim_, size_t horiz_dim_):  coord(new T[vert_dim_ * horiz_dim_]){
    transposed = false;
    vert_dim = vert_dim_;
    horiz_dim = horiz_dim_;
    if (v.size() != vert_dim_*horiz_dim_) {
        std::cerr<<"In Matrix(std::vector<T> v, size_t v, size_t h)  ";
        throw BadMatrixDimension();
    }
    for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
        coord[i] = v[i];
    }
}

template<class T>
Matrix<T>::Matrix(Matrix<T> &&rhs) noexcept {
    coord = rhs.coord;
    rhs.coord = nullptr;
    std::swap(transposed, rhs.transposed);
    std::swap(vert_dim, rhs.vert_dim);
    std::swap(horiz_dim, rhs.horiz_dim);
}

template<class T>
Matrix<T>::Matrix(const Matrix<T> &rhs) {
    transposed = rhs.transposed;
    vert_dim = rhs.vert_dim;
    horiz_dim = rhs.horiz_dim;
    coord = new T[rhs.vert_dim * rhs.horiz_dim];
    for (size_t i = 0; i < vert_dim * horiz_dim; i++) {
        coord[i] = rhs.coord[i];
    }
}

template<class T>
Matrix<T> &Matrix<T>::operator=(Matrix<T> &&rhs) noexcept {
    if (this == &rhs) return *this;
    Matrix<T> mat(std::move(rhs));
    std::swap(transposed, mat.transposed);
    std::swap(coord, mat.coord);
    std::swap(vert_dim, mat.vert_dim);
    std::swap(horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs) {
    if (this == &rhs) return *this;
    Matrix<T> mat(rhs);
    std::swap(transposed, mat.transposed);
    std::swap(coord, mat.coord);
    std::swap(vert_dim, mat.vert_dim);
    std::swap(horiz_dim, mat.horiz_dim);
    return *this;
}

template<class T>
Matrix<T>::~Matrix() {
    if (coord != nullptr)   {delete[] coord;}
}

template<class T>
std::vector<size_t> Matrix<T>::Shape() const {
    std::vector<size_t> shape({vert_dim, horiz_dim});
    return shape;
}

template<class T>
size_t Matrix<T>::HorizDim() const {
    return horiz_dim;
}

template<class T>
size_t Matrix<T>::VertDim() const {
    return vert_dim;
}

template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &rm) const {
    if (Shape() != rm.Shape())  {
        std::cerr<<"In Matrix::operator +   ";
        throw BadMatrixDimension();
    }
    T* res = new T[vert_dim*horiz_dim];
    for (size_t i = 0; i < vert_dim; ++i) {
        for (size_t j = 0; j < horiz_dim; ++j) {
            res[i*horiz_dim + j] = (*this)(i, j) + rm(i, j);
        }
    }
    return Matrix<T>(std::move(res), vert_dim, horiz_dim);
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &rm) const {
    if (Shape() != rm.Shape())  {
        std::cerr<<"In Matrix::operator -   ";
        throw BadMatrixDimension();
    }
    T* res = new T[vert_dim*horiz_dim];
    for (size_t i = 0; i < vert_dim; ++i) {
        for (size_t j = 0; j < horiz_dim; ++j) {
            res[i*horiz_dim + j] = (*this)(i, j) - rm(i, j);
        }
    }
    return Matrix<T>(std::move(res), vert_dim, horiz_dim);
}

template<class T>
Matrix<T> Matrix<T>::operator*(const T &val) const {
    size_t mat_size = vert_dim * horiz_dim;
    T* res = new T[mat_size];
    for (size_t i = 0; i < mat_size; i++) {
        res[i] = coord[i] * val;
    }
    return Matrix<T>(std::move(res), vert_dim, horiz_dim);
}

template<class T>
bool Matrix<T>::operator==(const Matrix<T> &rm) const {
    if (vert_dim != rm.vert_dim || horiz_dim != rm.horiz_dim) {
        return false;
    }
    for (size_t i = 0; i < vert_dim; ++i) {
        for (size_t j = 0; j < horiz_dim; j++) {
            if (std::abs((*this)(i, j) - rm(i, j)) > EPS) { return false; }
        }
    }
    return true;
}

template<class T>
bool Matrix<T>::operator!=(const Matrix<T> &rm) const {
    if (vert_dim != rm.vert_dim || horiz_dim != rm.horiz_dim) {
        return true;
    }
    for (size_t i = 0; i < vert_dim; ++i) {
        for (size_t j = 0; j < horiz_dim; j++) {
            if (std::abs((*this)(i, j) - rm(i, j)) > EPS) { return true; }
        }
    }
    return false;
}

template<class T>
const T &Matrix<T>::operator()(size_t i, size_t j) const {
    if (i < vert_dim && j < horiz_dim) {
        if (transposed)
            return coord[j * vert_dim + i];
        else
            return coord[i * horiz_dim + j];
    }
    else {
        throw IndexOutOfMatrix();
    }
}

template<class T>
T &Matrix<T>::operator()(size_t i, size_t j) {
    if (i < vert_dim && j < horiz_dim) {
        if (transposed)
            return coord[j * vert_dim + i];
        else
            return coord[i * horiz_dim + j];
    }
    else {
        throw IndexOutOfMatrix();
    }
}

template<class T>
Matrix<T> Matrix<T>::operator* (const Matrix<T>& rm) const {
    if (horiz_dim != rm.vert_dim) {
        std::cerr<<"In Matrix::operator *   ";
        throw BadMatrixDimension();
    }
    if (rm.horiz_dim * vert_dim * horiz_dim < 8000) {
        T *res = new T[rm.horiz_dim * vert_dim];
        T buffer;

        for (size_t i = 0; i < vert_dim; i++) {
            for (size_t j = 0; j < rm.horiz_dim; j++) {
                buffer = 0;
                for (size_t k = 0; k < horiz_dim; k++) {
                    buffer += (*this)(i, k) * rm(k, j);
                }
                res[i * rm.horiz_dim + j] = buffer;
            }
        }
        return Matrix(res, vert_dim, rm.horiz_dim);
    }
    else {
        Matrix<T> copied(rm);
        Matrix<T> this_copied(*this);
        return this_copied.SmartMult(copied);
    }
}

template<class T>
Matrix<T> Matrix<T>::SmartMult(Matrix<T> &rm) { //to fix non-8-divisible cases and transposed cases
    if (horiz_dim != rm.vert_dim) {
        std::cerr<<"In Matrix::SmartMult   ";
        throw BadMatrixDimension();
    }

    size_t M = (vert_dim/8 + 1) * 8;
    size_t K = (horiz_dim/8 + 1) * 8;
    size_t N = (rm.horiz_dim/8 + 1) * 8;

    T* left_mat = (*this).ExpandAndHardTranspose(M, K, false);
    T* right_mat = rm.ExpandAndHardTranspose(K, N, true);

    T* res = new T[M * N];
    T* res_norm = new T[vert_dim * rm.horiz_dim];
    gemm<T> multiply;
    multiply(M, K, N, left_mat, right_mat, res);

    for (size_t i = 0; i < vert_dim; ++i) {
        size_t cur_raw = i * rm.horiz_dim;
        size_t cur_8raw = i * N;
        for (size_t j = 0; j < rm.horiz_dim; ++j) {
            res_norm[cur_raw + j] = res[cur_8raw + j];
        }
    }
    return Matrix(res_norm, vert_dim, rm.horiz_dim);
}

template<class T>
T Matrix<T>::Determinant() const {
    if (VertDim() != HorizDim()) {
        std::cerr<<"In Matrix::Determinant   ";
        throw NotSquareMatrix();
    }
    T det_calc = T(1);
    std::vector<size_t> v;
    size_t rank = 0;
    Matrix<T> B = StraightRun(*this, Matrix<T>(horiz_dim, 1), det_calc, v, rank).first;

    for (size_t i = 0; i < B.VertDim(); i++) {
        det_calc *= B(i, i);
    }

    return det_calc;
}

template<class T>
size_t Matrix<T>::Rank() const {
    size_t rank = 0;
    T det_calc = T(0);
    std::vector<size_t> v;
    Matrix<T> B = StraightRun(*this, Matrix<T>(vert_dim, 1), det_calc, v, rank).first;
    return rank;
}

template<class T>
void Matrix<T>::Transpose() {
    transposed = !transposed;
    std::swap(vert_dim, horiz_dim);
}

template<class T>
T *Matrix<T>::ExpandAndHardTranspose(size_t M, size_t K, bool flip) {
    T* left_mat;
    left_mat = new T[M * K];
    if (!flip) {
        for (size_t i = 0; i < M; ++i) {
            size_t cur_raw = i * K;
            if (i < vert_dim) {
                for (size_t j = 0; j < K; ++j) {
                    if (j < horiz_dim)
                        left_mat[cur_raw + j] = (*this)(i, j);
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
                        left_mat[cur_col + i] = (*this)(i, j);
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

template<class T>
Matrix<T> Matrix<T>::Pow(Matrix<T> rm, int a) {
    int b=a;
    int i=0;
    while (b > 0) {
        b = b / 2;
        i++;
    }
    Matrix<T> a_2[i];//here we keep matrices for our degree
    int b_2[i];//here we keep numbers for our degree --4 is equal to [0,0,1]
    Matrix<T> degree[i];//here we keep matrices for the degrees of 2 [1,1,1,1,...]
    Matrix A(rm.VertDim(), rm.HorizDim(), IDENTITY);
    degree[0]=rm;
    //std::cout<<i<<" "<<b_2[0]<<" "<<b_2[1] <<b_2[2] <<std::endl;
    Matrix<T> Result(rm.VertDim(), rm.HorizDim(), IDENTITY);

    if (a % 2 == 1){a_2[0] = degree[0]*(a % 2);}
    else a_2[0] = A;
    a = a / 2;
    Result =a_2[0]*Result;
    if(i>1) {
        for (int k = 1; k < i; k++) {
            degree[k] = degree[k - 1] * degree[k - 1];
            if (a % 2 == 1) { a_2[k] = degree[k] * (a % 2); }
            else {
                a_2[k] = A;
            }
            a = a / 2;
            Result = a_2[k] * Result;
        }
    }
    return Result;
}

template<class T>
template<typename Q>
std::enable_if_t<is_complex<Q>::value, void> Matrix<T>::Conjugate() {
    for (size_t i = 0; i< VertDim(); i++) {
        for (size_t j = 0; j < HorizDim(); j++) {
            (*this)(i, j) = std::conj((*this)(i, j));
        }
    }
}

template<class T>
template<typename Q>
std::enable_if_t<!is_complex<Q>::value, void> Matrix<T>::Conjugate() {
    return;
}

template<typename T>
std::pair<bool, Matrix<T>> Inverse(const Matrix<T> & A) {
    if (A.VertDim() != A.HorizDim()) {
        std::cerr<<"In Inverse(const Matrix<T> &)   ";
        throw NotSquareMatrix();
    }
    T det_calc = T(1);
    size_t rank = 0;
    std::pair<Matrix<T>, Matrix<T>> p;
    std::vector<size_t> v;

    p = StraightRun(A, Matrix<T>(A.HorizDim(), A.HorizDim(), IDENTITY), det_calc, v, rank);

    if (det_calc == 0) {
        return {false, Matrix<T>(A.VertDim(), A.HorizDim(), NULLMATRIX)};
    }
    p = ReverseRun(p.first, p.second, v, rank);
    return {true, p.second};
}

template<typename T>
void ShowMatrix(const Matrix<T>&  mat)
{
    for (size_t i = 0; i < mat.VertDim(); i++) {
        for (size_t j = 0; j < mat.HorizDim(); j++) {
            printf(" %10.3f", mat(i,j));
        }
        printf("\n");
    }
    printf("\n");
}


#endif //LINALG_MATRIX_H