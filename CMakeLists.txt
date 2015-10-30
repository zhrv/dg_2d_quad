cmake_minimum_required(VERSION 3.0)
project(dg_2d_quad)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp)
add_executable(dg_2d_quad ${SOURCE_FILES})