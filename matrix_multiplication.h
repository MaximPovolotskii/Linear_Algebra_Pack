//
// Created by Максим on 24.11.2021.
//

#ifndef LINALG_MATRIX_MULTIPLICATION_H
#define LINALG_MATRIX_MULTIPLICATION_H
#include <immintrin.h>

template <typename U>
void gemm_v2(size_t M, size_t N, size_t K, const U* A, const U* B, U* C)
{
    for (size_t i = 0; i < M; ++i)
    {
        U* c = C + i * N;
        for (int j = 0; j < N; j += 8)
            _mm256_storeu_ps(c + j + 0, _mm256_setzero_ps());
        for (int k = 0; k < K; ++k)
        {
            const U *b = B + k * N;
            __m256 a = _mm256_set1_ps(A[i*K + k]);
            for (int j = 0; j < N; j += 16)
            {
                _mm256_storeu_ps(c + j + 0, _mm256_fmadd_ps(a,
                                                            _mm256_loadu_ps(b + j + 0), _mm256_loadu_ps(c + j + 0)));
                _mm256_storeu_ps(c + j + 8, _mm256_fmadd_ps(a,
                                                            _mm256_loadu_ps(b + j + 8), _mm256_loadu_ps(c + j + 8)));
            }
        }
    }
}

#endif //LINALG_MATRIX_MULTIPLICATION_H
