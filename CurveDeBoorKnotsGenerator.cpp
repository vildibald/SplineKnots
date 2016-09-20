#include "stdafx.h"
#include "CurveDeboorKnotsGenerator.h"
#include "SplineKnots.h"
#include "StopWatch.h"


splineknots::CurveDeboorKnotsGenerator::
CurveDeboorKnotsGenerator(const MathFunction function, bool optimized_tridiagonal)
	: function_(function), tridiagonal_(4,optimized_tridiagonal)
{
}

splineknots::CurveDeboorKnotsGenerator::
CurveDeboorKnotsGenerator(const InterpolativeMathFunction function, 
	bool optimized_tridiagonal)
	:function_(function), tridiagonal_(4, optimized_tridiagonal)
{
}

void splineknots::CurveDeboorKnotsGenerator::
InitializeKnots(const SurfaceDimension& dimension, KnotVector& knots)
{
	auto h = abs(dimension.max - dimension.min) / (dimension.knot_count - 1);
	auto x = dimension.min;
	for (auto i = 0; i < dimension.knot_count; i++ , x += h)
	{
		auto y = function_.Z()(x, 0);
		knots[i] = y;
	}
	tridiagonal_.ResizeBuffers(dimension.knot_count-2);
}

void splineknots::CurveDeboorKnotsGenerator::
RightSide(const KnotVector& knots, double h, double dfirst, double dlast)
{
	auto& rhs = tridiagonal_.RightSideBuffer();
	int num_unknowns = knots.size() - 2;
	auto mu1 = 3 / h;
	for (int i = 0; i < num_unknowns; i++)
	{
		rhs[i] = mu1 * (knots[i + 2] - knots[i]);
	}
	rhs[0] = rhs[0] - dfirst;
	rhs[num_unknowns - 1] = rhs[num_unknowns - 1] - dlast;
}

KnotVector splineknots::CurveDeboorKnotsGenerator::
GenerateKnots(const SurfaceDimension& dimension, double* calculation_time)
{
	StopWatch sw;
	
	KnotVector knots(dimension.knot_count);
	InitializeKnots(dimension, knots);
	auto dfirst = function_.Dx()(dimension.min, 0);
	auto dlast = function_.Dx()(dimension.max, 0);
	KnotVector result(knots.size());

	sw.Start();
	RightSide(knots, abs(dimension.max - dimension.min)
		/ (dimension.knot_count - 1), dfirst, dlast);
	
	auto& rhs = tridiagonal_.Solve(dimension.knot_count-2);
	result[0] = dfirst;
	result[result.size() - 1] = dlast;
	memcpy(&result.front() + 1, &rhs.front(), rhs.size());
	sw.Stop();
	
	if (calculation_time != nullptr)
	{
		*calculation_time = sw.EllapsedTime();
	}
	return result;
}

splineknots::Tridiagonal& splineknots::CurveDeboorKnotsGenerator::Tridiagonal()
{
	return tridiagonal_;
}