#pragma once
#include "stdafx.h"
#include <vector>
#include "Tridiagonal.h"


namespace splineknots
{
	class ReducedTridiagonal;
	typedef std::vector<ReducedTridiagonal> ReducedTridiagonals;
	class ReducedTridiagonal final
	{
		//friend class Tridiagonal;
		Tridiagonal tridiagonal_;
	public:

		ReducedTridiagonal(bool buffered = true);

		void ResizeBuffers(size_t newsize, bool shrinking_allowed = false);

		std::vector<double>& Solve(size_t num_unknowsns);
		
		std::vector<double>& RightSideBuffer();

	private:
		void ResizeBuffer(size_t newsize, bool shrinking_allowed = false);
		
		void ResizeRightSide(size_t newsize, bool shrinking_allowed = false);	
	public:
		std::vector<double>& ResetBufferAndGet();
		
		std::vector<double>& Buffer();
		
		const double& MainDiagonalValue() const;
	};
}