#include "stdafx.h"
#include "ReducedCurveDeboorGenerator.h"
#include "utils.h"
#include "StopWatch.h"


splineknots::ReducedCurveDeboorKnotsGenerator::ReducedCurveDeboorKnotsGenerator
(const MathFunction& function, bool optimized_tridiagonal)
	: function_(function), using_optimized_tridiagonal_(optimized_tridiagonal), 
	tridiagonal_(optimized_tridiagonal)
{
}

splineknots::ReducedCurveDeboorKnotsGenerator::ReducedCurveDeboorKnotsGenerator
(const InterpolativeMathFunction& function, bool optimized_tridiagonal) 
	: function_(function), using_optimized_tridiagonal_(optimized_tridiagonal),
	tridiagonal_(optimized_tridiagonal)
{
}

void splineknots::ReducedCurveDeboorKnotsGenerator::
RightSide(const KnotVector& function_values, double h, 
	double dfirst, double dlast)
{
	int n1, n2, N = function_values.size();
	n1 = (N % 2 == 0) ? (N - 2) / 2 : (N - 3) / 2;
	n2 = (N % 2 == 0) ? (N - 2) : (N - 3);

	auto& rhs = tridiagonal_.RightSideBuffer();
	double mu1 = 3.0 / h;
	double mu2 = 4.0 * mu1;
	int k = -1;
	for (int i = 2; i < n2; i+=2)
	{
		++k;
		rhs[k] = mu1 * (function_values[i + 1] - function_values[i - 2]) 
			- mu2 * (function_values[i + 1] - function_values[i - 1]);
	}
	rhs[0] -= dfirst;
	rhs[n1 - 1] -= dlast;
}

void splineknots::ReducedCurveDeboorKnotsGenerator::InitializeKnots(
	const SurfaceDimension& dimension, KnotVector& knots)
{
	auto h = abs(dimension.max - dimension.min) / (dimension.knot_count - 1);
	auto x = dimension.min;
	for (auto i = 0; i < dimension.knot_count; i++, x += h)
	{
		auto y = function_.Z()(x, 0);
		knots[i] = y;
	}
	tridiagonal_.ResizeBuffers(knots.size()/2 - 1);
}

KnotVector splineknots::ReducedCurveDeboorKnotsGenerator::
GenerateKnots(const splineknots::SurfaceDimension& dimension, 
	double* calculation_time)
{
	StopWatch sw;

	KnotVector knots(dimension.knot_count);
	InitializeKnots(dimension, knots);
	auto dfirst = function_.Dx()(dimension.min, 0);
	auto dlast = function_.Dx()(dimension.max, 0);
	auto h = abs(dimension.max - dimension.min) / (dimension.knot_count - 1);
	int n1, N = knots.size();
	n1 = (N % 2 == 0) ? (N - 2) / 2 : (N - 3) / 2;
	KnotVector result(knots.size());
	sw.Start();
	RightSide(knots, h, dfirst, dlast);
	auto rhs = tridiagonal_.Solve(n1);
	result[0] = dfirst;
	result[result.size() - 1] = dlast;

	for (int i = 0; i < n1; i++)
	{
		result[2 * (i + 1)] = rhs[i];
	}
	auto mu1 = 3.0 / h;
	for (int i = 1; i < result.size(); i += 2)
	{
		result[i] = 0.25 * (mu1 * (knots[i + 1] - knots[i - 1]) - 
			result[i - 1]
			- result[i + 1]);
	}
	sw.Stop();
	if (calculation_time != nullptr)
	{
		*calculation_time = sw.EllapsedTime();
	}
	return result;
}
