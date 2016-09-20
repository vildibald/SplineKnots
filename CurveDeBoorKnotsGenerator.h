#pragma once
#include <vector>
#include "MathFunction.h"
#include "SurfaceDimension.h"
#include "Tridiagonal.h"
#include "KnotVector.h"

namespace splineknots
{
	
	class CurveDeboorKnotsGenerator final
	{
		InterpolativeMathFunction function_;
		splineknots::Tridiagonal tridiagonal_;
		
		void RightSide(const KnotVector& function_values, double h, 
			double dfirst, double dlast);
		
		void InitializeKnots(const SurfaceDimension& dimension, 
			KnotVector& knots);
	public:
		CurveDeboorKnotsGenerator(const MathFunction function, 
			bool optimized_tridiagonal = true);

		CurveDeboorKnotsGenerator(const InterpolativeMathFunction function, 
			bool optimized_tridiagonal = true);

		KnotVector GenerateKnots(const SurfaceDimension& dimension, 
			double* calculation_time = nullptr);

		Tridiagonal& Tridiagonal();
	};
}
