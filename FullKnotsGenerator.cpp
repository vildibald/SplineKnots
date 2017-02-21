//
// Created by Viliam on 15.2.2017.
//
#include "stdafx.h"
#include "FullKnotsGenerator.h"
#include "utils.h"

splineknots::FullKnotsGenerator::FullKnotsGenerator(
        splineknots::MathFunction math_function, bool buffered)
        : function(math_function), tridagonals(), isParallel(false) {
    tridagonals.push_back(splineknots::Tridiagonal(4, buffered));
}

splineknots::FullKnotsGenerator::FullKnotsGenerator(
        splineknots::InterpolativeMathFunction math_function, bool buffered)
        : function(math_function),
          tridagonals(), isParallel(false) {
    tridagonals.push_back(splineknots::Tridiagonal(4, buffered));
}

void splineknots::FullKnotsGenerator::InParallel(bool value) {
    isParallel = value;
    auto threads = utils::num_threads;
    if (value) {
        tridagonals.reserve(threads);
        for (auto i = tridagonals.size(); i < threads; i++) {
            // create copy of tridiagonal solver
            auto copy_of_first = tridagonals[0];
            tridagonals.push_back(std::move(copy_of_first));
        }
    } else {
        for (int i = 0; i < tridagonals.size() - 1; ++i) {
            tridagonals.pop_back();
        }
    }
}

bool splineknots::FullKnotsGenerator::IsParallel() {
    return isParallel;
}

splineknots::Tridiagonals &splineknots::FullKnotsGenerator::Tridagonals() {
    return <#initializer#>;
}

Tridiagonal &splineknots::FullKnotsGenerator::Tridiagonal() {
    return <#initializer#>;
}

void splineknots::FullKnotsGenerator::InitializeBuffers(const size_t u_count,
                                                                  const size_t v_count) {

}

void splineknots::FullKnotsGenerator::Precalculate(
        const splineknots::SurfaceDimension &udimension,
        const splineknots::SurfaceDimension &vdimension) {

}

void splineknots::FullKnotsGenerator::InitializeKnots(
        const splineknots::SurfaceDimension &udimension,
        const splineknots::SurfaceDimension &vdimension, splineknots::KnotMatrix &values) {

}

void splineknots::FullKnotsGenerator::FillXDerivations(splineknots::KnotMatrix &values) {
    utils::For(0, static_cast<int>(values.ColumnsCount()),
               [&](int j) {
                   FillXDerivations(j, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillXYDerivations(splineknots::KnotMatrix &values) {
    FillXYDerivations(0, values);
    FillXYDerivations(values.ColumnsCount() - 1, values);
}

void splineknots::FullKnotsGenerator::FillYDerivations(splineknots::KnotMatrix &values) {
    utils::For(0, static_cast<int>(values.RowsCount()),
               [&](int i) {
                   FillYDerivations(i, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillYXDerivations(splineknots::KnotMatrix &values) {
    utils::For(0, static_cast<int>(values.RowsCount()),
               [&](int i) {
                   FillYXDerivations(i, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillXDerivations(const int column_index,
                                                                 splineknots::KnotMatrix &values) {
    auto unknowns_count = values.RowsCount() - 2;
    if (unknowns_count == 0) return;

    auto dset = [&](int index, double value) {
        values.SetDx(index, column_index, value);
    };
    auto rget = [&](int index) {
        return values.Z(index, column_index);
    };

    auto dlast = values.Dx(values.RowsCount() - 1, column_index);
    auto dfirst = values.Dx(0, column_index);

    SolveTridiagonal(rget, dfirst, dlast,
                     unknowns_count, dset);
}

void splineknots::FullKnotsGenerator::FillXYDerivations(const int column_index,
                                                                  splineknots::KnotMatrix &values) {

}

void splineknots::FullKnotsGenerator::FillYDerivations(const int row_index,
                                                                 splineknots::KnotMatrix &values) {

}

void splineknots::FullKnotsGenerator::FillYXDerivations(const int row_index,
                                                                  splineknots::KnotMatrix &values) {

}


