cmake_minimum_required(VERSION 3.6)
project(SplineKnots)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11 -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=core2")

set(SOURCE_FILES main.cpp
        CurveDeBoorKnotsGenerator.cpp CurveDeBoorKnotsGenerator.h
        ComparisonBenchmarkResult.cpp ComparisonBenchmarkResult.h
        FullKnotsGenerator.cpp FullKnotsGenerator.h
        EnhancedKnotsGenerator.cpp EnhancedKnotsGenerator.h
        KnotMatrix.cpp KnotMatrix.h
        KnotVector.h
        MathFunction.cpp MathFunction.h
        MulVsDiv.cpp MulVsDiv.h
        ReducedCurveDeboorGenerator.cpp ReducedCurveDeboorGenerator.h
        ReducedKnotsGenerator.cpp ReducedKnotsGenerator.h
        ReducedTridiagonal.cpp ReducedTridiagonal.h
        SplineKnots.h
        stdafx.cpp stdafx.h
        StopWatch.h
        SurfaceDimension.cpp SurfaceDimension.h
        Tridiagonal.cpp Tridiagonal.h
        utils.cpp utils.h utils_template.cpp
        NonUniformFullKnotsGenerator.cpp NonUniformFullKnotsGenerator.h)

add_executable(SplineKnots ${SOURCE_FILES})