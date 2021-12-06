//
// Created by Максим on 24.11.2021.
//

#ifndef LINALG_MATRIX_MULTIPLICATION_H
#define LINALG_MATRIX_MULTIPLICATION_H
#include <immintrin.h>

void float_gemm(size_t M, size_t N, size_t K, const float* A, const float* B, float* C)
{
    for (size_t i = 0; i < M; ++i)
    {
        float* c = C + i * N;
        for (int j = 0; j < N; j += 8)
            _mm256_storeu_ps(c + j + 0, _mm256_setzero_ps());
        for (int k = 0; k < K; ++k)
        {
            const float *b = B + k * N;
            __m256 a = _mm256_set1_ps(A[i*K + k]);
            for (int j = 0; j < N; j += 8)
            {
                _mm256_storeu_ps(c + j + 0, _mm256_fmadd_ps(a,
                                                            _mm256_loadu_ps(b + j + 0), _mm256_loadu_ps(c + j + 0)));
            }
        }
    }
}

void double_gemm(size_t M, size_t N, size_t K, const double* A, const double* B, double* C)
{
    for (size_t i = 0; i < M; ++i)
    {
        double* c = C + i * N;
        for (int j = 0; j < N; j += 8)
            _mm256_storeu_pd(c + j + 0, _mm256_setzero_pd());
        for (int k = 0; k < K; ++k)
        {
            const double *b = B + k * N;
            __m256d a = _mm256_set1_pd(A[i*K + k]);
            for (int j = 0; j < N; j += 8)
            {
                _mm256_storeu_pd(c + j + 0, _mm256_fmadd_pd(a,
                                                            _mm256_loadu_pd(b + j + 0), _mm256_loadu_pd(c + j + 0)));
            }
        }
    }
}

template <typename T, typename Enable = void>
struct gemm {
    void operator()(size_t M, size_t N, size_t K, const T* A, const T* B, T* C) {
        for (size_t i = 0; i < M; ++i) {
            size_t A_raw = i * K;
            for (size_t j = 0; j < N; ++j) {
                C[A_raw + j] = 0;
            }
            for (size_t k = 0; k < K; ++k) {
                size_t B_raw = k * M;
                size_t A_pos = A_raw + k;
                for (size_t j = 0; j < N; ++j) {
                    C[A_raw + j] += A[A_pos] * B[B_raw + j];
                }
            }
        }
    }
};

template<class T>
struct gemm <T, typename std::enable_if<std::is_same<T, float>::value>> {
    void operator()(size_t M, size_t N, size_t K, const T* A, const T* B, T* C) {
        float_gemm_v2(M, N, K, A, B, C);
    }
};

template<class T>
struct gemm <T, typename std::enable_if<std::is_same<T, double>::value>> {
    void operator()(size_t M, size_t N, size_t K, const T* A, const T* B, T* C) {
        double_gemm_v2(M, N, K, A, B, C);
    }
};

#endif //LINALG_MATRIX_MULTIPLICATION_H
