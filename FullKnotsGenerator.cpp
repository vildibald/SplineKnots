//
// Created by Viliam on 15.2.2017.
//
#include "stdafx.h"
#include "FullKnotsGenerator.h"
#include "utils.h"
#include "StopWatch.h"

splineknots::FullKnotsGenerator::FullKnotsGenerator(
        splineknots::MathFunction math_function, bool buffered)
        : function(math_function), uTridiagonals(), vTridiagonals(), isParallel(false) {

}

splineknots::FullKnotsGenerator::FullKnotsGenerator(
        splineknots::InterpolativeMathFunction math_function, bool buffered)
        : function(math_function),
          uTridiagonals(), vTridiagonals(), isParallel(false) {
}

void splineknots::FullKnotsGenerator::InParallel(bool value) {
    isParallel = value;
    auto threads = utils::logicalThreadCount;
    if (value) {
        uTridiagonals.reserve(threads);
        vTridiagonals.reserve(threads);
        for (auto i = uTridiagonals.size(); i < threads; i++) {
            // create copy of tridiagonal solver
            uTridiagonals.emplace_back(uTridiagonals[0]);
        }
        for (auto i = vTridiagonals.size(); i < threads; i++) {
            // create copy of tridiagonal solver
            vTridiagonals.emplace_back(vTridiagonals[0]);
        }
//    } else  {
//        for (int i = 0; i < uTridiagonals.size() - 1; ++i) {
//            uTridiagonals.pop_back();
//        }
//        for (int i = 0; i < vTridiagonals.size() - 1; ++i) {
//            vTridiagonals.pop_back();
//        }
    }
}


splineknots::KnotMatrix splineknots::FullKnotsGenerator::InitializeKnots(KnotVector uKnots,
                                                            KnotVector vKnots) {
    KnotMatrix values(std::move(uKnots), std::move(vKnots));

    // Init Z
    for (size_t i = 0; i < values.RowsCount(); i++) {
        for (size_t j = 0; j < values.ColumnsCount(); j++) {
            auto z = function.Z()(values.X(i), values.Y(j));
            //Function.Z(u,v); //Z(u, v);
            values.SetZ(i, j, z);
        }
    }
    // Init Dx
    auto uKnotCountMin1 = values.RowsCount() - 1;
    for (size_t j = 0; j < values.ColumnsCount(); j++) {
        auto dx = function.Dx()(values.X(j), values.Y(0));
        values.SetDx(0, j, dx); //Function.Dx(values(0,j).X, values(0,j).Y);
        values.SetDx(uKnotCountMin1, j, function.Dx()(values.X(j),
                                                      values.Y(uKnotCountMin1)));
    }
    // Init Dy
    auto vKnotCountMin1 = values.ColumnsCount() - 1;
    for (size_t i = 0; i < values.RowsCount(); i++) {
        values.SetDy(i, 0, function.Dy()(values.X(0), values.Y(i)));
        values.SetDy(i, vKnotCountMin1,
                     function.Dy()(values.X(vKnotCountMin1), values.Y(i))
        );
    }
    // Init Dxy
    values.SetDxy(0, 0, function.Dxy()(values.X(0), values.Y(0)));
    values.SetDxy(uKnotCountMin1, 0, function.Dxy()(values.X(0),
                                                    values.Y(uKnotCountMin1)));
    values.SetDxy(0, vKnotCountMin1, function.Dxy()(values.X(vKnotCountMin1),
                                                    values.Y(0)));
    values.SetDxy(uKnotCountMin1, vKnotCountMin1, function.Dxy()(values.X(vKnotCountMin1),
                                                                 values.Y(uKnotCountMin1)));
    return values;
}

void splineknots::FullKnotsGenerator::FillXDerivations(splineknots::KnotMatrix& values) {
    utils::For(0, static_cast<int>(values.ColumnsCount()),
               [&](int j) {
                   FillXDerivations(j, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillXYDerivations(splineknots::KnotMatrix& values) {
    FillXYDerivations(0, values);
    FillXYDerivations(values.ColumnsCount() - 1, values);
}

void splineknots::FullKnotsGenerator::FillYDerivations(splineknots::KnotMatrix& values) {
    utils::For(0, static_cast<int>(values.RowsCount()),
               [&](int i) {
                   FillYDerivations(i, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillYXDerivations(splineknots::KnotMatrix& values) {
    utils::For(0, static_cast<int>(values.RowsCount()),
               [&](int i) {
                   FillYXDerivations(i, values);
               },
               1, isParallel);
}

void splineknots::FullKnotsGenerator::FillXDerivations(const int column_index,
                                                       splineknots::KnotMatrix& values) {
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

    SolveTridiagonal(rget, dfirst, dlast, uTridiagonals[omp_get_thread_num()], dset);
}

void splineknots::FullKnotsGenerator::FillXYDerivations(const int column_index,
                                                        splineknots::KnotMatrix& values) {
    auto unknowns_count = values.RowsCount() - 2;
    if (unknowns_count == 0) return;

    auto dset = [&](int index, double value) {
        values.SetDxy(index, column_index, value);
    };
    auto rget = [&](int index) {
        return values.Dy(index, column_index);
    };

    auto dlast = values.Dxy(values.RowsCount() - 1, column_index);
    auto dfirst = values.Dxy(0, column_index);

    SolveTridiagonal(rget, dfirst, dlast, uTridiagonals[omp_get_thread_num()], dset);
}

void splineknots::FullKnotsGenerator::FillYDerivations(const int row_index,
                                                       splineknots::KnotMatrix& values) {
    auto unknowns_count = values.ColumnsCount() - 2;
    if (unknowns_count == 0) return;

    auto dset = [&](int index, double value) {
        values.SetDy(row_index, index, value);
    };
    auto rget = [&](int index) {
        return values.Z(row_index, index);
    };

    auto dlast = values.Dy(row_index, values.ColumnsCount() - 1);
    auto dfirst = values.Dy(row_index, 0);

    SolveTridiagonal(rget, dfirst, dlast, vTridiagonals[omp_get_thread_num()], dset);
}

void splineknots::FullKnotsGenerator::FillYXDerivations(const int row_index,
                                                        splineknots::KnotMatrix& values) {
    auto unknowns_count = values.ColumnsCount() - 2;
    if (unknowns_count == 0) return;

    auto dset = [&](int index, double value) {
        values.SetDxy(row_index, index, value);
    };
    auto rget = [&](int index) {
        return values.Dx(row_index, index);
    };

    auto dlast = values.Dxy(row_index, values.ColumnsCount() - 1);
    auto dfirst = values.Dxy(row_index, 0);

    SolveTridiagonal(rget, dfirst, dlast, vTridiagonals[omp_get_thread_num()], dset);
}

void splineknots::FullKnotsGenerator::InitializeTridiagonal(const KnotVector& uKnots,
                                                            const KnotVector& vKnots) {
    if(uTridiagonals.empty()) {
        uTridiagonals.emplace_back(uKnots, uKnots.size());
        vTridiagonals.emplace_back(vKnots, vKnots.size());
        return;
    }
    uTridiagonals[omp_get_thread_num()] = Tridiagonal(uKnots, uKnots.size());
    vTridiagonals[omp_get_thread_num()] = Tridiagonal(vKnots, uKnots.size());

}

splineknots::KnotMatrix
splineknots::FullKnotsGenerator::GenerateKnots(KnotVector uKnots,
                                               KnotVector vKnots,
                                               double* calculationTime) {
    StopWatch sw;

    if(uKnots.size()<4 || vKnots.size() < 4){
        return KnotMatrix::NullMatrix();
    }
    InitializeTridiagonal(uKnots, vKnots);
    auto values = InitializeKnots(std::move(uKnots), std::move(vKnots));

    sw.Start();
    FillXDerivations(values);
    FillXYDerivations(values);
    FillYDerivations(values);
    FillYXDerivations(values);
    sw.Stop();

    if (calculationTime != nullptr) {
        *calculationTime = sw.EllapsedTime();
    }
    return values;
}


