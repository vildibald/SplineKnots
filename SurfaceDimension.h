#pragma once
namespace splineknots
{
	struct SurfaceDimension
	{
		double min, max;
		int knotCount;
		double h;


		SurfaceDimension(double min, double max, int knot_count);

		double H() const;
	};
}
