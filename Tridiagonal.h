#pragma once

#include "stdafx.h"
#include <vector>
#include "KnotMatrix.h"
#include "KnotVector.h"


namespace splineknots {
    class Tridiagonal;

    typedef std::vector<Tridiagonal> Tridiagonals;

    class Tridiagonal final {
        std::vector<double> luBuffer;
        std::vector<double> rightSideBuffer;
        std::vector<double> deltas;
        std::vector<double> mainDiagonal;
        std::vector<double> deltaIMin1DivDeltaI;
        std::vector<double> deltaIDivDeltaIMin1;

        Tridiagonal(const std::vector<double>& deltas, const std::vector<double>& mainDiagonal,
                    const std::vector<double>& deltaIMin1DivDeltaI,
                    const std::vector<double>& deltaIDivDeltaIMin1,
                    const size_t numUnknowns
        );

    public:

        Tridiagonal(const KnotVector& knots, size_t numUnknowns);

        KnotVector& Solve();

        KnotVector& ResetBufferAndGet();

        KnotVector& Buffer();

        const std::vector<double>& GetDeltas() const;

        const std::vector<double>& GetDeltaIMin1DivDeltaI() const;

        const std::vector<double>& GetDeltaIDivDeltaIMin1() const;

        std::vector<double>& GetRightSideBuffer();

        //friends
//		friend class ReducedTridiagonal;
    };
}