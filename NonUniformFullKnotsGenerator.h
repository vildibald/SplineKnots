//
// Created by Viliam on 15.2.2017.
//
#include "MathFunction.h"
#include "KnotMatrix.h"
#include "SurfaceDimension.h"
#include <functional>
#include <vector>
#include "Tridiagonal.h"
#include <omp.h>


#ifndef SPLINEKNOTS_NONUNIFORMFULLKNOTSGENERATOR_H
#define SPLINEKNOTS_NONUNIFORMFULLKNOTSGENERATOR_H

namespace splineknots {
    class NonUniformFullKnotsGenerator {
        InterpolativeMathFunction function_;
        Tridiagonals tridagonals_;
        bool is_parallel_;
    public:

        NonUniformFullKnotsGenerator(MathFunction math_function, bool buffered = true);

        NonUniformFullKnotsGenerator(InterpolativeMathFunction math_function,
        bool buffered = true);

        KnotMatrix GenerateKnots(const SurfaceDimension& udimension,
                                 const SurfaceDimension& vdimension,
                                 double* calculation_time = nullptr);

        void InParallel(bool value);

        bool IsParallel();

        NonUniformFullKnotsGenerator(const MathFunction math_function,
                           Tridiagonal tridiagonal);

        NonUniformFullKnotsGenerator(const InterpolativeMathFunction math_function,
                           Tridiagonal tridiagonal);

        Tridiagonals& Tridagonals();

        Tridiagonal& Tridiagonal();

        void InitializeBuffers(const size_t u_count, const size_t v_count);

        void Precalculate(const SurfaceDimension& udimension,
                          const SurfaceDimension& vdimension);

        void InitializeKnots(const SurfaceDimension& udimension,
                             const SurfaceDimension& vdimension, KnotMatrix& values);

        void FillXDerivations(KnotMatrix& values);

        void FillXYDerivations(KnotMatrix& values);

        void FillYDerivations(KnotMatrix& values);

        void FillYXDerivations(KnotMatrix& values);

        void FillXDerivations(const int column_index, KnotMatrix& values);

        void FillXYDerivations(const int column_index, KnotMatrix& values);

        void FillYDerivations(const int row_index, KnotMatrix& values);

        void FillYXDerivations(const int row_index, KnotMatrix& values);

        template <typename RightSideSelector>
        void RightSide(const RightSideSelector& right_side_variables,
                       const double h, const double dfirst,
                       const double dlast, const int unknowns_count,
                       KnotVector& rightside)
        {
            auto h3 = 3/h;
            rightside[0] = h3 * (right_side_variables(2)
                                 - right_side_variables(0)) - dfirst;
            rightside[unknowns_count - 1] =
                    h3 * (right_side_variables(unknowns_count + 1)
                          - right_side_variables(unknowns_count - 1)) - dlast;
            for (auto i = 1; i < unknowns_count - 1; i++)
            {
                rightside[i] = h3 * (right_side_variables(i + 2)
                                     - right_side_variables(i));
            }
        }

        template <typename RightSideSelector, typename UnknownsSetter>
        void SolveTridiagonal(const RightSideSelector& selector,
                              const double h, const double dfirst,
                              const double dlast, const int unknowns_count,
                              UnknownsSetter& unknowns_setter)
        {
            auto& tridiagonal = Tridiagonal();
            auto& rightside = Tridiagonal().RightSideBuffer();
            RightSide(selector, h, dfirst, dlast, unknowns_count, rightside);

            auto& result = tridiagonal.Solve(unknowns_count);
            for (int k = 0; k < unknowns_count; k++)
            {
                unknowns_setter(k + 1, result[k]);
            }
        }
    };
}


#endif //SPLINEKNOTS_NONUNIFORMFULLKNOTSGENERATOR_H
