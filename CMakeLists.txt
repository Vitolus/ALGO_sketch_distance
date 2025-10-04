cmake_minimum_required(VERSION 3.22)
project(ALGO_sketch_distance)

set(CMAKE_CXX_STANDARD 17)

add_executable(dist main.cpp
        FracMinHash.cpp
        FracMinHash.h
        phylogenerator.h)
