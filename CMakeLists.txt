cmake_minimum_required(VERSION 3.17)
project(LAII)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-march=native")
add_executable(LAII catch_test.cpp qr_eigvals.h)