cmake_minimum_required(VERSION 3.22)
project(ALGO_sketch_distance)

set(CMAKE_CXX_STANDARD 17)

add_executable(ALGO_sketch_distance dist.cpp
        FracMinHash.cpp
        FracMinHash.h)
