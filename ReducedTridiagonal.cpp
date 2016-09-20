#include "stdafx.h"
#include <algorithm>
#include "SplineKnots.h"
#include "ReducedTridiagonal.h"
#include "Tridiagonal.h"


KnotVector&
splineknots::ReducedTridiagonal::Solve(size_t num_unknowns)
{
	auto even = num_unknowns % 2 == 0;
	auto num_equations = even ? num_unknowns / 2 - 1 : num_unknowns / 2;
	auto& rightside = RightSideBuffer();
	auto resize = std::max(num_equations, rightside.size());
	auto minsize = std::min(Buffer().size(), rightside.size());
	if (resize > minsize)
		ResizeBuffers(resize);
	double last_maindiag_value = even ? -15 : -14;
	auto& buffer = Buffer();
	utils::SolveDeboorTridiagonalSystemBuffered(MainDiagonalValue(),
		&rightside.front(), num_equations, &buffer.front(), 
		last_maindiag_value);
	return rightside;
}

splineknots::ReducedTridiagonal::ReducedTridiagonal(bool optimized)
	:tridiagonal_(-14, optimized)
{
}

void splineknots::ReducedTridiagonal::ResizeBuffers(size_t newsize,
                                                    bool shrinking_allowed)
{
	tridiagonal_.ResizeBuffers(newsize, shrinking_allowed);
}

KnotVector&
splineknots::ReducedTridiagonal::RightSideBuffer()
{
	return tridiagonal_.RightSideBuffer();
}

void splineknots::ReducedTridiagonal::ResizeBuffer(size_t newsize,
                                                   bool shrinking_allowed)
{
	tridiagonal_.ResizeBuffer(newsize, shrinking_allowed);
}

void splineknots::ReducedTridiagonal::ResizeRightSide(size_t newsize,
                                                      bool shrinking_allowed)
{
	tridiagonal_.ResizeRightSide(newsize, shrinking_allowed);
}

KnotVector&
splineknots::ReducedTridiagonal::ResetBufferAndGet()
{
	return tridiagonal_.ResetBufferAndGet();
}

KnotVector&
splineknots::ReducedTridiagonal::Buffer()
{
	return tridiagonal_.Buffer();
}


const double& splineknots::ReducedTridiagonal::MainDiagonalValue() const
{
	return tridiagonal_.MainDiagonalValue();
}
