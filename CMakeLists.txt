cmake_minimum_required(VERSION 3.6)
project(SplineKnots)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++14 -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=core2")

set(SOURCE_FILES main.cpp
        ComparisonBenchmarkResult.cpp ComparisonBenchmarkResult.h
        Spline.cpp Spline.h
        MathFunction.cpp MathFunction.h
        OperationsBenchmark.cpp OperationsBenchmark.h
        stdafx.cpp stdafx.h
        StopWatch.h
        SurfaceDimension.cpp SurfaceDimension.h
        utils.cpp utils.h
        SplineKnots.h FullAlgorithm.cpp FullAlgorithm.h ReducedAlgorithm.cpp ReducedAlgorithm.h)

add_executable(SplineKnots ${SOURCE_FILES})