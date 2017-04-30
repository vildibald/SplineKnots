#pragma once
#include "stdafx.h"
#include <vector>
#include "KnotMatrix.h"


namespace splineknots
{
	class Tridiagonal;
	typedef std::vector<Tridiagonal> Tridiagonals;
	class Tridiagonal final
	{	
		std::vector<double> lu_buffer_;
		std::vector<double> right_side_buffer_;
		const double main_diagonal_value;
	public:
		
		Tridiagonal(double main_value, bool optimized = true);

		void ResizeBuffers(size_t newsize, bool shrinking_allowed = false);

		std::vector<double>& Solve(size_t num_unknowsns);
		std::vector<double>& RightSideBuffer();
		void ResizeBuffer(size_t newsize, bool shrinking_allowed = false);
		void ResizeRightSide(size_t newsize, bool shrinking_allowed = false);
		std::vector<double>& ResetBufferAndGet();
		std::vector<double>& Buffer();
		const double& MainDiagonalValue() const;
		bool IsUsingOptimizedTridiagonal() const;

		//friends
		friend class ReducedTridiagonal;
	};
}