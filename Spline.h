#pragma once

#include <vector>
#include "SurfaceDimension.h"

namespace splineknots
{
	class Spline
	{
		SurfaceDimension xdim;
		SurfaceDimension ydim;
		std::vector<std::vector<double>> z;
		std::vector<std::vector<double>> dx;
		std::vector<std::vector<double>> dy;
		std::vector<std::vector<double>> dxy;

		Spline();
	public:
		static Spline NullMatrix();
		
		bool IsNull();

		Spline(SurfaceDimension rowdimension, SurfaceDimension columndimension);
		
		size_t RowsCount() const
		{
			return xdim.knotCount;
		}

		size_t ColumnsCount() const
		{
			return ydim.knotCount;
		}

		double X(const size_t j) const
		{
			auto h = abs(xdim.max - xdim.min) / (xdim.knotCount-1);
			return j*h + xdim.min;
		}

		double Y(const size_t i) const
		{
			auto h = abs(ydim.max - ydim.min) / (ydim.knotCount-1);
			return i*h + ydim.min;
		}

		const SurfaceDimension& XDimensionParameters() const
		{
			return xdim;
		}

		const SurfaceDimension& YDimensionParameters() const
		{
			return ydim;
		}

		double Z(const size_t i, const size_t j) const
		{
			return z[i][j];
		}
		double Dx(const size_t i, const size_t j) const
		{
			return dx[j][i];
		}
		double Dy(const size_t i, const size_t j) const
		{
			return dy[i][j];
		}
		double Dxy(const size_t i, const size_t j) const
		{
			return dxy[i][j];
		}
	
		void SetZ(const size_t i, const size_t j, const double value)
		{
			z[i][j] = value;
		}
		void SetDx(const size_t i, const size_t j, const double value)
		{
			dx[j][i] = value;
		}
		void SetDy(const size_t i, const size_t j, const double value)
		{
			dy[i][j] = value;
		}
		void SetDxy(const size_t i, const size_t j, const double value)
		{
			dxy[i][j] = value;
		}

		void Print();
	};
}
