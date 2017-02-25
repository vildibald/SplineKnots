cmake_minimum_required(VERSION 3.6)
project(SplineKnots)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=core2")

set(SOURCE_FILES main.cpp
        ComparisonBenchmarkResult.cpp ComparisonBenchmarkResult.h
        FullKnotsGenerator.cpp FullKnotsGenerator.h
        KnotMatrix.cpp KnotMatrix.h
        KnotVector.h
        MathFunction.cpp MathFunction.h
        MulVsDiv.cpp MulVsDiv.h
        SplineKnots.h
        stdafx.cpp stdafx.h
        StopWatch.h
        Tridiagonal.cpp Tridiagonal.h
        utils.cpp utils.h
        )

add_executable(SplineKnots ${SOURCE_FILES})