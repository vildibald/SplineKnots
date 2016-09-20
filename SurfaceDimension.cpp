#include "stdafx.h"
#include "SurfaceDimension.h"


namespace splineknots
{
	SurfaceDimension::SurfaceDimension(double min, double max, int knotCount)
		:min(min),
		 max(max),
		 knot_count(knotCount)
	{
	}
}
