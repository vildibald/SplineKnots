#pragma once
#include "CurveDeboorKnotsGenerator.h"
#include "ReducedTridiagonal.h"
namespace splineknots
{
	class ReducedCurveDeboorKnotsGenerator final
	{
		InterpolativeMathFunction function_;
		bool using_optimized_tridiagonal_;
		ReducedTridiagonal tridiagonal_;

		void RightSide(const KnotVector& function_values, 
			double h, double dfirst, double dlast);
		
		void InitializeKnots(const SurfaceDimension& dimension,
			KnotVector& knots);
	public:
		KnotVector GenerateKnots(const SurfaceDimension& dimension, 
			double* calculation_time = nullptr);
		
		ReducedCurveDeboorKnotsGenerator(const MathFunction& function, 
			bool optimized_tridiagonal = true);
		
		ReducedCurveDeboorKnotsGenerator
			(const InterpolativeMathFunction& function, bool optimized_tridiagonal = true);
	};
}
