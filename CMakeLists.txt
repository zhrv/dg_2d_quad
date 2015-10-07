cmake_minimum_required(VERSION 3.3)
project(dg_2d_quad)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp GaussIntegrator1d.cpp GaussIntegrator2D.cpp)
add_executable(dg_2d_quad ${SOURCE_FILES})