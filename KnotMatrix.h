#pragma once

#include "KnotVector.h"

namespace splineknots
{
	class KnotMatrix
	{

		double** z;
		double** dx;
		double** dy;
		double** dxy;
		KnotVector x;
		KnotVector y;

		KnotMatrix();
	public:
		static KnotMatrix NullMatrix();
		
		bool IsNull();
		
		KnotMatrix(KnotVector rowVector, KnotVector columnVector);
		
		KnotMatrix(const KnotMatrix& other);
		
		KnotMatrix(KnotMatrix&& other) noexcept;
		
		KnotMatrix& operator =(const KnotMatrix& other);
		
		KnotMatrix& operator =(KnotMatrix&& other) noexcept;
		
		~KnotMatrix() noexcept;

		size_t RowsCount() const
		{
			return x.size();
		}

		size_t ColumnsCount() const
		{
			return y.size();
		}


		double** Z() const
		{
			return z;
		}

		double** Dx() const
		{
			return dx;
		}

		double** Dy() const
		{
			return dy;
		}

		double** Dxy() const
		{
			return dxy;
		}

		const KnotVector& X() const{
			return x;
		}

		const KnotVector& Y() const{
			return y;
		}

		const double X(size_t i) const{
			return x[i];
		}

		const double Y(size_t j) const{
			return y[j];
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
