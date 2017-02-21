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


#ifndef SPLINEKNOTS_FULLKNOTSGENERATOR_H
#define SPLINEKNOTS_FULLKNOTSGENERATOR_H

namespace splineknots {
    class FullKnotsGenerator {
        InterpolativeMathFunction function;
        Tridiagonals tridiagonals;
        bool isParallel;
    public:

        FullKnotsGenerator(MathFunction mathFunction, bool buffered = true);

        FullKnotsGenerator(InterpolativeMathFunction mathFunction,
                                     bool buffered = true);

        KnotMatrix GenerateKnots(const SurfaceDimension &udimension,
                                 const SurfaceDimension &vdimension,
                                 double *calculationTime = nullptr);

        void InParallel(bool value);

        bool IsParallel();

        FullKnotsGenerator(const MathFunction mathFunction,
                                     Tridiagonal tridiagonal);

        FullKnotsGenerator(const InterpolativeMathFunction mathFunction,
                                     Tridiagonal tridiagonal);

        Tridiagonals &Tridagonals();

        Tridiagonal &Tridiagonal();

        void InitializeBuffers(const size_t u_count, const size_t v_count);

        void Precalculate(const SurfaceDimension &udimension,
                          const SurfaceDimension &vdimension);

        void InitializeKnots(const SurfaceDimension &udimension,
                             const SurfaceDimension &vdimension, KnotMatrix &values);

        void FillXDerivations(KnotMatrix &values);

        void FillXYDerivations(KnotMatrix &values);

        void FillYDerivations(KnotMatrix &values);

        void FillYXDerivations(KnotMatrix &values);

        void FillXDerivations(const int column_index, KnotMatrix &values);

        void FillXYDerivations(const int column_index, KnotMatrix &values);

        void FillYDerivations(const int row_index, KnotMatrix &values);

        void FillYXDerivations(const int row_index, KnotMatrix &values);

        template<typename RightSideSelector, typename KnotCoordinateSelector>
        void RightSide(const RightSideSelector &rightSideSelector,
                       const KnotCoordinateSelector &knotCoordinateSelector,
                       const double dfirst,
                       const double dlast, const int unknownsCount,
                       KnotVector &rightSide) {

            auto deltaXIPlus1 = knotCoordinateSelector(2) - knotCoordinateSelector(1);
            auto deltaXI = knotCoordinateSelector(1) - knotCoordinateSelector(0);
            rightSide[0] = 3 * (
                    deltaXI / deltaXIPlus1 * (rightSideSelector(2) - rightSideSelector(1))
                    +
                    deltaXIPlus1 / deltaXI * (rightSideSelector(1) - rightSideSelector(0))
            );
            deltaXIPlus1 = knotCoordinateSelector(2) - knotCoordinateSelector(1);
            deltaXI = knotCoordinateSelector(1) - knotCoordinateSelector(0);
            rightSide[unknownsCount - 1] = 3 * (
                    deltaXI / deltaXIPlus1 * (rightSideSelector(unknownsCount + 1)
                                              - rightSideSelector(unknownsCount))
                    +
                    deltaXIPlus1 / deltaXI * (rightSideSelector(unknownsCount)
                                              - rightSideSelector(unknownsCount - 1))
            );
            for (auto i = 1; i < unknownsCount - 1; i++) {
                deltaXIPlus1 = knotCoordinateSelector(i + 2) - knotCoordinateSelector(i + 1);
                deltaXI = knotCoordinateSelector(i + 1) - knotCoordinateSelector(i);
                rightSide[i] = 3 * (
                        deltaXI / deltaXIPlus1 * (rightSideSelector(i + 2) - rightSideSelector
                                (i + 1))
                        +
                        deltaXIPlus1 / deltaXI * (rightSideSelector(i + 1) - rightSideSelector(i))
                );
            }
        }

        template<typename RightSideSelector, typename KnotCoordinateSelector,
                typename UnknownsSetter>
        void SolveTridiagonal(const RightSideSelector &selector,
                              const KnotCoordinateSelector& knotCoordinateSelector,
                              const double dfirst,
                              const double dlast, const int unknowns_count,
                              UnknownsSetter &unknowns_setter) {
            auto &tridiagonal = Tridiagonal();
            auto &rightside = Tridiagonal().RightSideBuffer();
            RightSide(selector, knotCoordinateSelector, dfirst, dlast, unknowns_count, rightside);

            auto &result = tridiagonal.Solve();
            for (int k = 0; k < unknowns_count; k++) {
                unknowns_setter(k + 1, result[k]);
            }
        }
    };
}


#endif //SPLINEKNOTS_FULLKNOTSGENERATOR_H