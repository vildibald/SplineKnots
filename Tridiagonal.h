#pragma once
#include "stdafx.h"
#include <vector>
#include "KnotMatrix.h"
#include "KnotVector.h"


namespace splineknots
{
	class Tridiagonal;
	typedef std::vector<Tridiagonal> Tridiagonals;
	class Tridiagonal final
	{	
		KnotVector lu_buffer_;
		KnotVector right_side_buffer_;
		const double main_diagonal_value;
	public:
		
		Tridiagonal(double main_value, bool optimized = true);

		void ResizeBuffers(size_t newsize, bool shrinking_allowed = false);

		KnotVector& Solve(size_t num_unknowsns);
		KnotVector& RightSideBuffer();
		void ResizeBuffer(size_t newsize, bool shrinking_allowed = false);
		void ResizeRightSide(size_t newsize, bool shrinking_allowed = false);
		KnotVector& ResetBufferAndGet();
		KnotVector& Buffer();
		const double& MainDiagonalValue() const;
		bool IsUsingOptimizedTridiagonal() const;

		//friends
		friend class ReducedTridiagonal;
	};
}