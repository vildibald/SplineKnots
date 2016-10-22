#include "stdafx.h"
#include "EnhancedReducedDeboorKnotsGenerator.h"
#include "StopWatch.h"
#include "utils.h"

splineknots::EnhancedReducedDeboorKnotsGenerator::EnhancedReducedDeboorKnotsGenerator(
        MathFunction math_function, bool buffered)
        : function_(math_function), tridagonals_(), is_parallel_(false),
          deboor_(math_function),
          precalculated_hx_(1), precalculated_hy_(1) {
    tridagonals_.push_back(ReducedTridiagonal(buffered));
}

splineknots::EnhancedReducedDeboorKnotsGenerator::EnhancedReducedDeboorKnotsGenerator(
        InterpolativeMathFunction math_function, bool buffered)
        : function_(math_function), tridagonals_(), is_parallel_(false),
          deboor_(math_function),
          precalculated_hx_(1), precalculated_hy_(1) {
    tridagonals_.push_back(ReducedTridiagonal(buffered));
}


void splineknots::EnhancedReducedDeboorKnotsGenerator::RightSideCross(
        const KnotMatrix &knots, const int i, const double dfirst,
        const double dlast, const int unknowns_count, KnotVector &rightside) {

}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillXDerivations(KnotMatrix &
values) {
    auto h = precalculated_hx_.deboor_precalculated_.h;
    auto three_div_h = precalculated_hx_.deboor_precalculated_.three_div_h;
    int nrows = values.RowsCount();
    int ncols = values.ColumnsCount();
    utils::For(0, ncols,
               [&](int j) {
                   FillXDerivations(j, values);
                   for (size_t i = 1; i < nrows - 1; i += 2) {
                       values.SetDx(i, j,
                                    0.25 * (three_div_h * (values.Z(i + 1, j) -
                                                           values.Z(i - 1, j)) -
                                            values.Dx(i + 1, j) +
                                            values.Dx(i - 1, j))
                       );
                   }
               },
               1, is_parallel_);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillXYDerivations(KnotMatrix &
values) {
    auto h = precalculated_hx_.deboor_precalculated_.h;
    auto three_div_h = precalculated_hx_.deboor_precalculated_.three_div_h;
    int nrows = values.RowsCount();
    int ncols = values.ColumnsCount();
    FillXYDerivations(0, values);
    FillXYDerivations(ncols - 1, values);
    for (size_t i = 1; i < nrows - 1; i += 2) {
        values.SetDxy(i, 0,
                      0.25 *
                      (three_div_h * (values.Dx(i + 1, 0) -
                                      values.Dx(i - 1, 0)) -
                       values.Dxy(i + 1, 0) +
                       values.Dxy(i - 1, 0))
        );
        values.SetDxy(i, ncols - 1,
                      0.25 *
                      (three_div_h * (values.Dx(i + 1, ncols - 1) -
                                      values.Dx(i - 1, ncols - 1)) -
                       values.Dxy(i + 1, ncols - 1) +
                       values.Dxy(i - 1, ncols - 1))
        );
    }


}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillYDerivations(KnotMatrix &
values) {
    auto h = precalculated_hy_.deboor_precalculated_.h;
    auto three_div_h = precalculated_hy_.deboor_precalculated_.three_div_h;
    int ncols = values.ColumnsCount();
    int nrows = values.RowsCount();
    utils::For(0, nrows,
               [&](int i) {
                   FillYDerivations(i, values);
                   for (size_t j = 1; j < ncols - 1; j += 2) {
                       values.SetDy(i, j,
                                    0.25 * (three_div_h * (values.Z(i, j + 1) -
                                                           values.Z(i, j - 1)) -
                                            values.Dy(i, j + 1) +
                                            values.Dy(i, j - 1))
                       );
                   }
               }, 1, is_parallel_);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillYXDerivations(KnotMatrix &
values) {
    auto h = precalculated_hy_.deboor_precalculated_.h;
    auto three_div_h = precalculated_hy_.deboor_precalculated_.three_div_h;
    int ncols = values.ColumnsCount();
    int nrows = values.RowsCount();
    utils::For(0, nrows,
               [&](int i) {
                   FillYXDerivations(i, values);
                   for (size_t j = 1; j < ncols - 1; j += 2) {
                       values.SetDxy(i, j,
                                     0.25 *
                                     (three_div_h * (values.Dy(i, j + 1) -
                                                     values.Dy(i, j - 1)) -
                                      values.Dxy(i, j + 1) +
                                      values.Dxy(i, j - 1))
                       );
                   }
               }, 1, is_parallel_);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillXDerivations(const int
                                                                   column_index,
                                                                   KnotMatrix &values) {
    auto unknowns_count = values.RowsCount() - 2;
    if (unknowns_count == 0)
        return;

    auto dset = [&](int index, double value) {
        values.SetDx(index, column_index, value);
    };
    auto rget = [&](int index) {
        return values.Z(index, column_index);
    };

    auto dlast = values.Dx(values.RowsCount() - 1, column_index);
    auto dfirst = values.Dx(0, column_index);

    SolveTridiagonal(rget, precalculated_hx_, dfirst, dlast, unknowns_count,
                     dset);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillXYDerivations(const int
                                                                    column_index,
                                                                    KnotMatrix &values) {
    auto unknowns_count = values.RowsCount() - 2;
    if (unknowns_count == 0)
        return;

    auto dset = [&](int index, double value) {
        values.SetDxy(index, column_index, value);
    };
    auto rget = [&](int index) {
        return values.Dy(index, column_index);
    };

    auto dlast = values.Dxy(values.RowsCount() - 1, column_index);
    auto dfirst = values.Dxy(0, column_index);

    SolveTridiagonal(rget, precalculated_hx_, dfirst, dlast, unknowns_count,
                     dset);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillYDerivations(const int
                                                                   row_index,
                                                                   KnotMatrix &values) {
    auto unknowns_count = values.ColumnsCount() - 2;
    if (unknowns_count == 0)
        return;

    auto dset = [&](int index, double value) {
        values.SetDy(row_index, index, value);
    };
    auto rget = [&](int index) {
        return values.Z(row_index, index);
    };

    auto dlast = values.Dy(row_index, values.ColumnsCount() - 1);
    auto dfirst = values.Dy(row_index, 0);

    SolveTridiagonal(rget, precalculated_hy_, dfirst, dlast, unknowns_count,
                     dset);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::FillYXDerivations(const int
                                                                    row_index,
                                                                    KnotMatrix &values) {
    auto unknowns_count = values.ColumnsCount() - 2;
    if (unknowns_count == 0)
        return;

    auto dset = [&](int index, double value) {
        values.SetDxy(row_index, index, value);
    };
    auto rget = [&](int index) {
        return values.Dy(row_index, index);
    };

    auto dlast = values.Dxy(row_index, values.ColumnsCount() - 1);
    auto dfirst = values.Dxy(row_index, 0);

    SolveTridiagonal(rget, precalculated_hy_, dfirst, dlast, unknowns_count,
                     dset);
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::InitializeBuffers(const size_t
                                                                    u_count,
                                                                    const size_t v_count) {
    auto size = std::max(u_count / 2 - 1, v_count / 2 - 1);
    auto &trid = tridagonals_;
    for (size_t i = 0; i < trid.size(); i++) {
        trid[i].ResizeBuffers(size);
    }
    deboor_.InitializeBuffers(u_count, v_count);
}


splineknots::KnotMatrix splineknots::EnhancedReducedDeboorKnotsGenerator::
GenerateKnots(const SurfaceDimension &udimension, const SurfaceDimension &
vdimension, double *calculation_time) {
    StopWatch sw;

    if (udimension.knot_count < 6 || vdimension.knot_count < 6) {
        deboor_.InitializeBuffers(udimension.knot_count,
                                  vdimension.knot_count);
        return deboor_.GenerateKnots(udimension, vdimension);
    }
    KnotMatrix values{udimension, vdimension};
    InitializeBuffers(udimension.knot_count, vdimension.knot_count);
    InitializeKnots(udimension, vdimension, values);

    sw.Start();
    FillXDerivations(values);
    FillYDerivations(values);
    FillXYDerivations(values);
    FillYXDerivations(values);
    sw.Stop();

    if (calculation_time != nullptr) {
        *calculation_time = sw.EllapsedTime();
    }
    return values;
}

void
splineknots::EnhancedReducedDeboorKnotsGenerator::InParallel(bool in_parallel) {
    is_parallel_ = in_parallel;
    auto threads = utils::num_threads;
    if (in_parallel) {
        tridagonals_.reserve(threads);
        for (auto i = tridagonals_.size(); i < threads; i++) {
            // create copy of tridiagonal solver
            auto copy_of_first = tridagonals_[0];
            tridagonals_.push_back(std::move(copy_of_first));
        }
    } else {
        for (int i = 0; i < tridagonals_.size() - 1; ++i) {
            tridagonals_.pop_back();
        }
    }
}

splineknots::ReducedTridiagonals &
splineknots::EnhancedReducedDeboorKnotsGenerator::Tridagonals() {
    return tridagonals_;
}

splineknots::ReducedTridiagonal &
splineknots::EnhancedReducedDeboorKnotsGenerator::Tridiagonal() {
    return tridagonals_[omp_get_thread_num()];
}


void splineknots::EnhancedReducedDeboorKnotsGenerator::InitializeKnots(
        const SurfaceDimension &udimension, const SurfaceDimension &vdimension,
        KnotMatrix &values) {
    Precalculate(udimension, vdimension);
    deboor_.InitializeKnots(udimension, vdimension, values);
}

void splineknots::EnhancedReducedDeboorKnotsGenerator::Precalculate(
        const SurfaceDimension &udimension,
        const SurfaceDimension &vdimension) {
    precalculated_hx_ = PrecalculatedEnhancedReduced(abs(
            udimension.max - udimension.min) / (udimension.knot_count - 1));
    precalculated_hy_ = PrecalculatedEnhancedReduced(abs(
            vdimension.max - vdimension.min) / (vdimension.knot_count - 1));
}

splineknots::EnhancedReducedDeboorKnotsGenerator::PrecalculatedEnhancedReduced
::PrecalculatedEnhancedReduced(const double h) : deboor_precalculated_(h) {
    twelve_div_h = 4 * deboor_precalculated_.three_div_h;
    three_h_div_4 = 0.75 * h;
}
