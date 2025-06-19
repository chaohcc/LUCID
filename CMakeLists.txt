cmake_minimum_required(VERSION 3.21)
project(lucid_model_guided_buffer)
set(CMAKE_CXX_STANDARD 17)


include(FetchContent)
include(ExternalProject)

#find_package(OpenMP)
#find_package(MKL REQUIRED)
#find_package(JeMalloc REQUIRED)

#include_directories(${MKL_INCLUDE_DIRS})

#find_package(TBB REQUIRED)
#include_directories(${JEMALLOC_INCLUDE_DIR})

#/usr/local/lib/gcc/x86_64-pc-linux-gnu/7.2.0/include/popcntintrin.h:42:1: error：inlining failed in call to always_inline ‘long long int _mm_popcnt_u64(long long unsigned int)’: target specific option mismatch
#_mm_popcnt_u64 (unsigned long long __X)   # this error can be solved by

if(MSVC)
    set(CMAKE_CXX_FLAGS "/O2 /arch:AVX2 /W1 /EHsc")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "-O3 -xHost")
else()
    # clang and gcc
    set(CMAKE_CXX_FLAGS "-O0 -march=native -Wall -Wextra")    #chao changes O3 to O0
endif()



add_executable(lucid_model_guided_buffer main.cpp lucid.h utils.h workload.h)
