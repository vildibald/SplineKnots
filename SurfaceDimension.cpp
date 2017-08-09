#include "stdafx.h"
#include "SurfaceDimension.h"


namespace splineknots
{
	SurfaceDimension::SurfaceDimension(double min, double max, int knotCount)
		:min(min),
		 max(max),
		 knotCount(knotCount),
		 h(abs(max-min)/(knotCount-1))
	{
	}

	double SurfaceDimension::H() const {
		return h;
	}
}
