#pragma once
namespace splineknots
{
	struct SurfaceDimension
	{
		double min, max;
		int knot_count;
		SurfaceDimension(double min, double max, int knot_count);

		SurfaceDimension(const SurfaceDimension& other) = default;
		SurfaceDimension(SurfaceDimension&& other) = default;
		SurfaceDimension& operator=(const SurfaceDimension& other) = default;
		SurfaceDimension& operator=(SurfaceDimension&& other) = default;
	};
}
