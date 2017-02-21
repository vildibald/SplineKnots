#pragma once

#include "stdafx.h"
#include <vector>
#include "KnotMatrix.h"
#include "KnotVector.h"


namespace splineknots {
    class Tridiagonal;

    typedef std::vector<Tridiagonal> Tridiagonals;

    class Tridiagonal final {
        KnotVector luBuffer;
        KnotVector rightSideBuffer;
        std::vector<double> lowerDiagonal_;
        std::vector<double> mainDiagonal_;
        std::vector<double> upperDiagonal_;

        Tridiagonal(std::vector<double> lowerDiagonal,
                    std::vector<double> mainDiagonal,
                    std::vector<double> upperDiagonal);

    public:
        template<typename KnotCoordinates>
        static Tridiagonal
        Create(KnotCoordinates &knots, size_t numUnknowns) {
            std::vector<double> lowerDiagonal(numUnknowns);
            std::vector<double> mainDiagonal(numUnknowns);
            std::vector<double> upperDiagonal(numUnknowns);

            for (int i = 1; i < numUnknowns - 1; ++i) {
                lowerDiagonal.emplace_back(
                        knots(i) - knots(i - 1)
                );
                mainDiagonal.emplace_back(2 * (
                        knots(i + 1) - knots(i - 1)
                ));
                upperDiagonal.emplace_back(
                        knots(i + 1) - knots(i)
                );
            }

            return Tridiagonal(lowerDiagonal, mainDiagonal, upperDiagonal);
        }


        void ResizeBuffers(size_t newsize, bool shrinking_allowed = false);

        KnotVector &Solve();

        KnotVector &RightSideBuffer();

        void ResizeBuffer(size_t newsize, bool shrinking_allowed = false);

        void ResizeRightSide(size_t newsize, bool shrinking_allowed = false);

        KnotVector &ResetBufferAndGet();

        KnotVector &Buffer();

        bool IsUsingOptimizedTridiagonal() const;

        //friends
//		friend class ReducedTridiagonal;
    };
}