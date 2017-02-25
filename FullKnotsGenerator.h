//
// Created by Viliam on 15.2.2017.
//
#include "MathFunction.h"
#include "KnotMatrix.h"
#include <functional>
#include <vector>
#include "Tridiagonal.h"
#include <omp.h>


#ifndef SPLINEKNOTS_FULLKNOTSGENERATOR_H
#define SPLINEKNOTS_FULLKNOTSGENERATOR_H

namespace splineknots {
    class FullKnotsGenerator {
        InterpolativeMathFunction function;
        Tridiagonals uTridiagonals;
        Tridiagonals vTridiagonals;
        bool isParallel;

        void InitializeTridiagonal(const KnotVector& uKnots,
                                   const KnotVector& vKnots);

        KnotMatrix InitializeKnots(KnotVector uKnots,
                                   KnotVector vKnots);

    public:

        FullKnotsGenerator(MathFunction mathFunction,
                           bool buffered = true);

        FullKnotsGenerator(InterpolativeMathFunction mathFunction,
                           bool buffered = true);

        KnotMatrix GenerateKnots(KnotVector uKnots,
                                 KnotVector vKnots,
                                 double* calculationTime = nullptr);

        void InParallel(bool value);

        bool IsParallel();


        void FillXDerivations(KnotMatrix& values);

        void FillXYDerivations(KnotMatrix& values);

        void FillYDerivations(KnotMatrix& values);

        void FillYXDerivations(KnotMatrix& values);

        void FillXDerivations(const int column_index,
                              KnotMatrix& values);

        void FillXYDerivations(const int column_index,
                               KnotMatrix& values);

        void FillYDerivations(const int row_index,
                              KnotMatrix& values);

        void FillYXDerivations(const int row_index,
                               KnotMatrix& values);

        template<typename RightSideSelector>
        void RightSide(const RightSideSelector& rightSideSelector,
                       const double dfirst,
                       const double dlast,
                       Tridiagonal& tridiagonal) {

            auto& rightSide = tridiagonal.GetRightSideBuffer();
            auto& deltaIMin1DivDeltaI = tridiagonal.GetDeltaIMin1DivDeltaI();
            auto& deltaIDivDeltaIMin1 = tridiagonal.GetDeltaIDivDeltaIMin1();
            auto& deltas = tridiagonal.GetDeltas();
            rightSide[0] = 3 * (
                    deltaIMin1DivDeltaI[0] * (rightSideSelector(2) - rightSideSelector(1))
                    +
                    deltaIDivDeltaIMin1[0] * (rightSideSelector(1) - rightSideSelector(0))
            ) - deltas[0] * dfirst;
            rightSide.back() = 3 * (
                    deltaIMin1DivDeltaI.back() * (rightSideSelector(rightSide.size()+1)
                                                  - rightSideSelector(rightSide.size()))
                    +
                    deltaIDivDeltaIMin1.back() * (rightSideSelector(rightSide.size())
                                                  - rightSideSelector(rightSide.size() - 1))
            ) - deltas.back() * dlast;
            for (auto i = 1; i < rightSide.size() - 1; i++) {
                rightSide[i] = 3 * (
                        deltaIMin1DivDeltaI[i] * (rightSideSelector(i + 2) - rightSideSelector
                                (i + 1))
                        +
                        deltaIDivDeltaIMin1[i] * (rightSideSelector(i + 1) - rightSideSelector(i))
                );
            }
        }

        template<typename RightSideSelector, typename UnknownsSetter>
        void SolveTridiagonal(const RightSideSelector& selector,
                              const double dfirst,
                              const double dlast,
                              Tridiagonal& tridiagonal,
                              UnknownsSetter& unknowns_setter) {
            RightSide(selector, dfirst, dlast, tridiagonal);

            auto& result = tridiagonal.Solve();
            for (int k = 0; k < result.size(); k++) {
                unknowns_setter(k + 1, result[k]);
            }
        }
    };
}


#endif //SPLINEKNOTS_FULLKNOTSGENERATOR_H