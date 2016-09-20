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

		KnotVector& Solve(size_t num_unknowsns);
		
		KnotVector& RightSideBuffer();

	private:
		void ResizeBuffer(size_t newsize, bool shrinking_allowed = false);
		
		void ResizeRightSide(size_t newsize, bool shrinking_allowed = false);	
	public:
		KnotVector& ResetBufferAndGet();
		
		KnotVector& Buffer();
		
		const double& MainDiagonalValue() const;
	};
}