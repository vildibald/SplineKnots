#pragma once
#include "SurfaceDimension.h"

namespace splineknots
{
	class KnotMatrix
	{
		SurfaceDimension xdim_;
		SurfaceDimension ydim_;
		double** z_;
		double** dx_;
		double** dy_;
		double** dxy_;

		KnotMatrix();
	public:
		static KnotMatrix NullMatrix();
		
		bool IsNull();
		
		KnotMatrix(SurfaceDimension rowdimension, SurfaceDimension columndimension);
		
		KnotMatrix(const KnotMatrix& other);
		
		KnotMatrix(KnotMatrix&& other) noexcept;
		
		KnotMatrix& operator =(const KnotMatrix& other);
		
		KnotMatrix& operator =(KnotMatrix&& other) noexcept;
		
		~KnotMatrix() noexcept;

		size_t RowsCount() const
		{
			return xdim_.knot_count;
		}

		size_t ColumnsCount() const
		{
			return ydim_.knot_count;
		}


		double** Z() const
		{
			return z_;
		}

		double** Dx() const
		{
			return dx_;
		}

		double** Dy() const
		{
			return dy_;
		}

		double** Dxy() const
		{
			return dxy_;
		}

		double X(const size_t i,const size_t j) const
		{
			return j*abs(xdim_.max - xdim_.min) / xdim_.knot_count;
		}

		double Y(const size_t i, const size_t j) const
		{
			return i*abs(ydim_.max - ydim_.min) / ydim_.knot_count;
		}

		const SurfaceDimension& XDimensionParameters() const
		{
			return xdim_;
		}

		const SurfaceDimension& YDimensionParameters() const
		{
			return ydim_;
		}

		double Z(const size_t i, const size_t j) const
		{
			return z_[i][j];
		}
		double Dx(const size_t i, const size_t j) const
		{
			return dx_[j][i];
		}
		double Dy(const size_t i, const size_t j) const
		{
			return dy_[i][j];
		}
		double Dxy(const size_t i, const size_t j) const
		{
			return dxy_[i][j];
		}
	
		void SetZ(const size_t i, const size_t j, const double value)
		{
			z_[i][j] = value;
		}
		void SetDx(const size_t i, const size_t j, const double value)
		{
			dx_[j][i] = value;
		}
		void SetDy(const size_t i, const size_t j, const double value)
		{
			dy_[i][j] = value;
		}
		void SetDxy(const size_t i, const size_t j, const double value)
		{
			dxy_[i][j] = value;
		}

		void Print();
	};
}
