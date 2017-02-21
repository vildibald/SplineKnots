//
// Created by Viliam on 15.2.2017.
//
#include "stdafx.h"
#include "NonUniformFullKnotsGenerator.h"

splineknots::NonUniformFullKnotsGenerator::NonUniformFullKnotsGenerator(
        splineknots::MathFunction math_function, bool buffered)
        : function_(math_function), tridagonals_(), is_parallel_(false) {
    tridagonals_.push_back(splineknots::Tridiagonal(4, buffered));
}

splineknots::NonUniformFullKnotsGenerator::NonUniformFullKnotsGenerator(
        splineknots::InterpolativeMathFunction math_function, bool buffered)
        : function_(math_function),
          tridagonals_(), is_parallel_(false) {
    tridagonals_.push_back(splineknots::Tridiagonal(4, buffered));
}



