cmake_minimum_required(VERSION 3.22)
project(ALGO_sketch_distance)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the OpenMP package
find_package(OpenMP REQUIRED)

add_executable(dist main.cpp
        FracMinHash.cpp
        FracMinHash.h
        phylogenerator.h
        BloomFilter.cpp
        BloomFilter.h)

# Add optimization flags
target_compile_options(dist PRIVATE -O3)

# Link against OpenMP. CMake will automatically add the required compiler flags.
if(OpenMP_FOUND)
    target_link_libraries(dist PUBLIC OpenMP::OpenMP_CXX)
endif()