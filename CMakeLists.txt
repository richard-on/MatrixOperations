cmake_minimum_required(VERSION 3.16.3)
project(mv)

set(CMAKE_CXX_STANDARD 20)

add_executable(mv main.cpp matrix.cpp matrix.h vector.cpp vector.h lu.cpp lu.h ldlt.cpp ldlt.h)
