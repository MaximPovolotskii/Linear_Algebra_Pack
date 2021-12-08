cmake_minimum_required(VERSION 3.17)
project(linalg)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-march=native")
add_executable(linalg main.cpp  matrix.h vector.h gauss_run.h linear_system_solve.h qr_decomposition.h)
add_executable(hh new_main.cpp matrix.h vector.h gauss_run.h linear_system_solve.h qr_decomposition.h)
add_executable(ei eigenvectors_test.cpp matrix.h vector.h gauss_run.h linear_system_solve.h qr_decomposition.h)
