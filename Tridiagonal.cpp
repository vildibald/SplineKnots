#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>


KnotVector&
splineknots::Tridiagonal::ResetBufferAndGet() {
    auto& buffer = luBuffer;
    std::fill(buffer.begin(), buffer.end(), 1);
    return buffer;
}

KnotVector&
splineknots::Tridiagonal::Buffer() {
    return luBuffer;
}

KnotVector&
splineknots::Tridiagonal::Solve() {
    auto numUnknowns = mainDiagonal.size();
    auto& buffer = Buffer();
    utils::SolveTridiagonalSystemBuffered(&deltas.front(), &mainDiagonal.front(),
                                          (&deltas.front()) + 1,
                                          &rightSideBuffer.front(), numUnknowns,
                                          &buffer.front());
    return rightSideBuffer;
}

splineknots::Tridiagonal::Tridiagonal(const std::vector<double>& deltas,
                                      const std::vector<double>& mainDiagonal,
                                      const std::vector<double>& deltaIMin1DivDeltaI,
                                      const std::vector<double>& deltaIDivDeltaIMin1,
                                      const size_t numUnknowns)
        : deltas(deltas),
          mainDiagonal(mainDiagonal),
          deltaIMin1DivDeltaI(deltaIMin1DivDeltaI),
          deltaIDivDeltaIMin1(deltaIDivDeltaIMin1),
          luBuffer(numUnknowns),
          rightSideBuffer(numUnknowns) {
    luBuffer.assign(numUnknowns, 0);
    rightSideBuffer.assign(numUnknowns, 0);
}

const std::vector<double>& splineknots::Tridiagonal::GetDeltas() const {
    return deltas;
}

const std::vector<double>& splineknots::Tridiagonal::GetDeltaIMin1DivDeltaI() const {
    return deltaIMin1DivDeltaI;
}

const std::vector<double>& splineknots::Tridiagonal::GetDeltaIDivDeltaIMin1() const {
    return deltaIDivDeltaIMin1;
}

std::vector<double>& splineknots::Tridiagonal::GetRightSideBuffer() {
    return rightSideBuffer;
}

splineknots::Tridiagonal::Tridiagonal(const KnotVector& knots, size_t numUnknowns) {
    std::vector<double> deltas;
    std::vector<double> mainDiagonal;
    std::vector<double> deltaIMin1DivDeltaI;
    std::vector<double> deltaIDivDeltaIMin1;
    numUnknowns-=2;
    deltas.reserve(numUnknowns);
    mainDiagonal.reserve(numUnknowns);
    deltaIMin1DivDeltaI.reserve(numUnknowns);
    deltaIDivDeltaIMin1.reserve(numUnknowns);

    for (int i = 0; i < numUnknowns; ++i) {
        deltas.emplace_back(
                knots[i + 1] - knots[i]
        );
    }

    for (int i = 1; i < numUnknowns+1; ++i) {
        mainDiagonal.emplace_back(2 * (
                deltas[i - 1] + deltas[i]
        ));

        deltaIMin1DivDeltaI.emplace_back(
                deltas[i - 1] / deltas[i]
        );

        deltaIDivDeltaIMin1.emplace_back(
                1 / deltaIMin1DivDeltaI.back()
        );
    }
    this->deltas = std::move(deltas);
    this->mainDiagonal = std::move(mainDiagonal);
    this->deltaIMin1DivDeltaI = std::move(deltaIMin1DivDeltaI);
    this->deltaIDivDeltaIMin1 = std::move(deltaIDivDeltaIMin1);
    luBuffer = std::vector<double>(numUnknowns);
    rightSideBuffer = std::vector<double>(numUnknowns);
}

