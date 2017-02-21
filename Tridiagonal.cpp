#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>


void splineknots::Tridiagonal::ResizeBuffers(size_t newsize, bool
shrinking_allowed) {
    ResizeBuffer(newsize, shrinking_allowed);
    ResizeRightSide(newsize, shrinking_allowed);
}

void splineknots::Tridiagonal::ResizeBuffer(size_t newsize,
                                            bool shrinking_allowed) {
    auto oldsize = luBuffer.size();
    if (newsize > oldsize || shrinking_allowed) {
        luBuffer.resize(newsize);
    }
}

void splineknots::Tridiagonal::ResizeRightSide(size_t newsize, bool
shrinking_allowed) {
    auto oldsize = rightSideBuffer.size();
    if (newsize > oldsize || shrinking_allowed) {
        rightSideBuffer.resize(newsize);
    }
}


KnotVector &
splineknots::Tridiagonal::ResetBufferAndGet() {
    auto &buffer = luBuffer;
    std::fill(buffer.begin(), buffer.end(), 1);
    return buffer;
}

KnotVector &
splineknots::Tridiagonal::Buffer() {
    return luBuffer;
}

KnotVector &
splineknots::Tridiagonal::Solve() {
    auto numUnknowns = mainDiagonal_.size();
    auto resize = std::max(numUnknowns, rightSideBuffer.size());
    auto minsize = std::min(Buffer().size(), RightSideBuffer().size());
    if (resize > minsize)
        ResizeBuffers(resize);
    auto &buffer = Buffer();
    utils::SolveTridiagonalSystemBuffered(&lowerDiagonal_.front(), &mainDiagonal_.front(),
                                          &upperDiagonal_.front(),
                                          &rightSideBuffer.front(), numUnknowns,
                                          &buffer.front());
    return rightSideBuffer;
}

KnotVector &
splineknots::Tridiagonal::RightSideBuffer() {
    return rightSideBuffer;
}


splineknots::Tridiagonal::Tridiagonal(std::vector<double> lowerDiagonal,
                                      std::vector<double> mainDiagonal,
                                      std::vector<double> upperDiagonal)
        : lowerDiagonal_(lowerDiagonal_), mainDiagonal_(mainDiagonal_),
          upperDiagonal_(upperDiagonal_) {

}


