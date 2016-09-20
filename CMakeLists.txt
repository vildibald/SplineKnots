cmake_minimum_required(VERSION 3.6)
project(SplineKnots)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=core2")

set(SOURCE_FILES main.cpp)
add_executable(SplineKnots ${SOURCE_FILES})