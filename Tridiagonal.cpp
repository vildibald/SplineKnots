#include "stdafx.h"
#include "Tridiagonal.h"
#include "utils.h"
#include <algorithm>


splineknots::Tridiagonal::Tridiagonal(double main_value,
                                      bool buffered)
	:lu_buffer_(),
	 right_side_buffer_(),
	 main_diagonal_value(main_value)
{
}

void splineknots::Tridiagonal::ResizeBuffers(size_t newsize, bool
                                             shrinking_allowed)
{
	ResizeBuffer(newsize, shrinking_allowed);
	ResizeRightSide(newsize, shrinking_allowed);
}

void splineknots::Tridiagonal::ResizeBuffer(size_t newsize,
                                            bool shrinking_allowed)
{
	auto oldsize = lu_buffer_.size();
	if (newsize > oldsize || shrinking_allowed)
	{
		lu_buffer_.resize(newsize);
	}
}

void splineknots::Tridiagonal::ResizeRightSide(size_t newsize, bool
                                               shrinking_allowed)
{
	auto oldsize = right_side_buffer_.size();
	if (newsize > oldsize || shrinking_allowed)
	{
		right_side_buffer_.resize(newsize);
	}
}


KnotVector&
splineknots::Tridiagonal::ResetBufferAndGet()
{
	auto& buffer = lu_buffer_;
	std::fill(buffer.begin(), buffer.end(), 1);
	return buffer;
}

KnotVector&
splineknots::Tridiagonal::Buffer()
{
	return lu_buffer_;
}

KnotVector&
splineknots::Tridiagonal::Solve(size_t num_unknowns)
{
	auto resize = std::max(num_unknowns, right_side_buffer_.size());
	auto minsize = std::min(Buffer().size(), RightSideBuffer().size());
	if (resize > minsize)
		ResizeBuffers(resize);
	auto& buffer = Buffer();
	utils::SolveDeboorTridiagonalSystemBuffered(
		main_diagonal_value,
		&right_side_buffer_.front(), num_unknowns, &buffer.front());
	return right_side_buffer_;
}

KnotVector&
splineknots::Tridiagonal::RightSideBuffer()
{
	return right_side_buffer_;
}

const double& splineknots::Tridiagonal::MainDiagonalValue() const
{
	return main_diagonal_value;
}
