//#include "stdafx.h"
//#include "EnhancedReducedDeboorKnotsGenerator.h"
//#include "StopWatch.h"
//
//splineknots::EnhancedReducedDeboorKnotsGenerator::
//EnhancedReducedDeboorKnotsGenerator(MathFunction function, bool buffered)
//	: function_(function), tridagonals_(), is_parallel_(false),
//	deboor_(function),
//	precalculated_hx_(1), precalculated_hy_(1)
//{
//}
//
//splineknots::EnhancedReducedDeboorKnotsGenerator::
//EnhancedReducedDeboorKnotsGenerator(InterpolativeMathFunction function, bool buffered)
//	: function_(function), tridagonals_(), is_parallel_(false),
//	deboor_(function),
//	precalculated_hx_(1), precalculated_hy_(1)
//{
//	tridagonals_.push_back(splineknots::Tridiagonal(4,buffered));
//}
//
//splineknots::KnotMatrix splineknots::EnhancedReducedDeboorKnotsGenerator::
//GenerateKnots(const SurfaceDimension& udimension,
//	const SurfaceDimension& vdimension, double* calculation_time)
//{
//	StopWatch sw;
//
//	if (udimension.knot_count < 6 || vdimension.knot_count < 6)
//	{
//		deboor_.InitializeBuffers(udimension.knot_count,
//			vdimension.knot_count);
//		return deboor_.GenerateKnots(udimension, vdimension);
//	}
//	KnotMatrix values{ udimension, vdimension };
//	InitializeBuffers(udimension.knot_count, vdimension.knot_count);
//	InitializeKnots(udimension, vdimension, values);
//
//	sw.Start();
//	FillXDerivations(values);
//	FillYDerivations(values);
//	FillXYDerivations(values);
//	FillYXDerivations(values);
//	sw.Stop();
//
//	if (calculation_time != nullptr)
//	{
//		*calculation_time = sw.EllapsedTime();
//	}
//	return values;
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::InParallel(bool in_parallel)
//{
//}
//
//splineknots::Tridiagonals& splineknots::EnhancedReducedDeboorKnotsGenerator::Tridagonals()
//{
//}
//
//splineknots::Tridiagonal& splineknots::EnhancedReducedDeboorKnotsGenerator::Tridiagonal()
//{
//}
//
//splineknots::EnhancedReducedDeboorKnotsGenerator::PrecalculatedEnhancedReduced::PrecalculatedEnhancedReduced(const double h)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::InitializeKnots(const SurfaceDimension& udimension, const SurfaceDimension& vdimension, KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::Precalculate(const SurfaceDimension& udimension, const SurfaceDimension& vdimension)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::InitializeBuffers(size_t u_count, size_t v_count)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::RightSideCross(const KnotMatrix& knots, const int i, const double dfirst, const double dlast, const int unknowns_count, KnotVector& rightside)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillXDerivations(KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillXYDerivations(KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillYDerivations(KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillYXDerivations(KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillXDerivations(const int column_index, KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillXYDerivations(const int column_index, KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillYDerivations(const int row_index, KnotMatrix& values)
//{
//}
//
//void splineknots::EnhancedReducedDeboorKnotsGenerator::FillYXDerivations(const int row_index, KnotMatrix& values)
//{
//}
