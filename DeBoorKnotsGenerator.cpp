#include "stdafx.h"
#include "DeBoorKnotsGenerator.h"

#include "Tridiagonal.h"

#include "SplineKnots.h"

#include "StopWatch.h"

namespace splineknots
{
	Tridiagonals& DeBoorKnotsGenerator::Tridagonals()
	{
		return tridagonals_;
	}

	splineknots::Tridiagonal& DeBoorKnotsGenerator::Tridiagonal()
	{
		return  tridagonals_[omp_get_thread_num()];
	}

	DeBoorKnotsGenerator::Precalculated::Precalculated(const double h)
		: h(h), three_div_h(3 / h)
	{
	}

	DeBoorKnotsGenerator::
		DeBoorKnotsGenerator(MathFunction math_function, bool buffered) :
		function_(math_function), tridagonals_(), is_parallel_(false), 
		precalculated_hx_(1), precalculated_hy_(1)
	{
		tridagonals_.push_back(splineknots::Tridiagonal(4, 
			buffered));
	}

	DeBoorKnotsGenerator::DeBoorKnotsGenerator(
		InterpolativeMathFunction math_function, bool buffered)
		: function_(math_function),
		tridagonals_(), is_parallel_(false),
		precalculated_hx_(1), precalculated_hy_(1)
	{
		tridagonals_.push_back(splineknots::Tridiagonal(4, buffered));
	}


	DeBoorKnotsGenerator::DeBoorKnotsGenerator(MathFunction math_function,
		splineknots::Tridiagonal tridiagonal)
		: function_(math_function),
		tridagonals_(), is_parallel_(false),
		precalculated_hx_(1), precalculated_hy_(1)
	{
		tridagonals_.push_back(std::move(tridiagonal));
	}

	DeBoorKnotsGenerator::DeBoorKnotsGenerator(InterpolativeMathFunction
		math_function, splineknots::Tridiagonal tridiagonal)
		: function_(math_function),
		tridagonals_(), is_parallel_(false),
		precalculated_hx_(1), precalculated_hy_(1)
	{
		tridagonals_.push_back(std::move(tridiagonal));
	}

	KnotMatrix DeBoorKnotsGenerator::GenerateKnots(const SurfaceDimension&
		udimension, const SurfaceDimension& vdimension, double*
		calculation_time)
	{
		StopWatch sw;
		
		if (udimension.knot_count < 4 || vdimension.knot_count < 4)
		{
			return KnotMatrix::NullMatrix();
		}
		KnotMatrix values{udimension, vdimension};
		InitializeBuffers(udimension.knot_count, vdimension.knot_count);
		InitializeKnots(udimension, vdimension, values);

		sw.Start();
		FillXDerivations(values);
		FillXYDerivations(values);
		FillYDerivations(values);
		FillYXDerivations(values);
		sw.Stop();

		if (calculation_time != nullptr)
		{
			*calculation_time = sw.EllapsedTime();
		}
		return values;
	}

	void DeBoorKnotsGenerator::Precalculate(const SurfaceDimension& udimension,
		const SurfaceDimension& vdimension)
	{
		SetPrecalculatedHX(Precalculated(abs(udimension.max
			- udimension.min) / (udimension.knot_count - 1)));
		SetPrecalculatedHY(Precalculated(abs(vdimension.max
			- vdimension.min) / (vdimension.knot_count - 1)));
	}

	void DeBoorKnotsGenerator::InitializeKnots(const SurfaceDimension&
		udimension, const SurfaceDimension& vdimension, KnotMatrix& values)
	{
		Precalculate(udimension, vdimension);
		auto hx = precalculated_hx_.h;
		auto hy = precalculated_hy_.h;
		auto u = udimension.min;
		auto& f = function_;
		// Init Z
		for (auto i = 0; i < udimension.knot_count; i++, u += hx)
		{
			auto v = vdimension.min;
			for (auto j = 0; j < vdimension.knot_count; j++, v += hy)
			{
				auto z = f.Z()(u, v);
				//Function.Z(u,v); //Z(u, v);
				values.SetZ(i, j, z);
			}
		}
		// Init Dx
		auto uKnotCountMin1 = udimension.knot_count - 1;
		for (auto j = 0; j < vdimension.knot_count; j++)
		{
			auto dx = f.Dx()(values.X(0, j), values.Y(0, j));
			values.SetDx(0, j, dx); //Function.Dx(values(0,j).X, values(0,j).Y);
			values.SetDx(uKnotCountMin1, j, f.Dx()(values.X(uKnotCountMin1, j),
				values.Y(uKnotCountMin1, j)));
		}
		// Init Dy
		auto vKnotCountMin1 = vdimension.knot_count - 1;
		for (auto i = 0; i < udimension.knot_count; i++)
		{
			values.SetDy(i, 0, f.Dy()(values.X(i, 0), values.Y(i, 0)));
			values.SetDy(i, vKnotCountMin1,
				f.Dy()(values.X(i, vKnotCountMin1), values.Y(i,
					vKnotCountMin1))
				);
		}
		// Init Dxy
		values.SetDxy(0, 0, f.Dxy()(values.X(0, 0), values.Y(0, 0)));
		values.SetDxy(uKnotCountMin1, 0, f.Dxy()(values.X(uKnotCountMin1, 0),
			values.Y(uKnotCountMin1, 0)));
		values.SetDxy(0, vKnotCountMin1, f.Dxy()(values.X(0, vKnotCountMin1),
			values.Y(0, vKnotCountMin1)));
		values.SetDxy(uKnotCountMin1, vKnotCountMin1, f.Dxy()(values.X(
			uKnotCountMin1, vKnotCountMin1), values.Y(uKnotCountMin1,
				vKnotCountMin1)));
	}

	void DeBoorKnotsGenerator::FillXDerivations(KnotMatrix& values)
	{
		utils::For(0, static_cast<int>(values.ColumnsCount()),
			[&](int j)
		{
			FillXDerivations(j, values);
		},
			1, is_parallel_);
	}

	void DeBoorKnotsGenerator::FillXYDerivations(KnotMatrix& values)
	{
		FillXYDerivations(0, values);
		FillXYDerivations(values.ColumnsCount() - 1, values);
	}

	void DeBoorKnotsGenerator::FillYDerivations(KnotMatrix& values)
	{
		utils::For(0, static_cast<int>(values.RowsCount()),
			[&](int i)
		{
			FillYDerivations(i, values);
		},
			1, is_parallel_);
	}

	void DeBoorKnotsGenerator::FillYXDerivations(KnotMatrix& values)
	{
		utils::For(0, static_cast<int>(values.RowsCount()),
			[&](int i)
		{
			FillYXDerivations(i, values);
		},
			1, is_parallel_);
	}

	void DeBoorKnotsGenerator::FillXDerivations(const int column_index,
		KnotMatrix& values)
	{
		auto unknowns_count = values.RowsCount() - 2;
		if (unknowns_count == 0) return;

		auto dset = [&](int index, double value)
		{
			values.SetDx(index, column_index, value);
		};
		auto rget = [&](int index)
		{
			return values.Z(index, column_index);
		};

		auto dlast = values.Dx(values.RowsCount() - 1, column_index);
		auto dfirst = values.Dx(0, column_index);

		SolveTridiagonal(rget, PrecalculatedHX(), dfirst, dlast,
			unknowns_count, dset);
	}

	void DeBoorKnotsGenerator::FillXYDerivations(const int column_index,
		KnotMatrix& values)
	{
		auto unknowns_count = values.RowsCount() - 2;
		if (unknowns_count == 0) return;

		auto dset = [&](int index, double value)
		{
			values.SetDxy(index, column_index, value);
		};
		auto rget = [&](int index)
		{
			return values.Dy(index, column_index);
		};

		auto h = precalculated_hx_.h;
		auto dlast = values.Dxy(values.RowsCount() - 1, column_index);
		auto dfirst = values.Dxy(0, column_index);

		SolveTridiagonal(rget, PrecalculatedHX(), dfirst, dlast,
			unknowns_count, dset);
	}

	void DeBoorKnotsGenerator::FillYDerivations(const int row_index,
		KnotMatrix& values)
	{
		auto unknowns_count = values.ColumnsCount() - 2;
		if (unknowns_count == 0) return;

		auto dset = [&](int index, double value)
		{
			values.SetDy(row_index, index, value);
		};
		auto rget = [&](int index)
		{
			return values.Z(row_index, index);
		};

		auto dlast = values.Dy(row_index, values.ColumnsCount() - 1);
		auto dfirst = values.Dy(row_index, 0);

		SolveTridiagonal(rget, PrecalculatedHY(), dfirst, dlast,
			unknowns_count, dset);
	}

	void DeBoorKnotsGenerator::FillYXDerivations(const int row_index,
		KnotMatrix& values)
	{
		auto unknowns_count = values.ColumnsCount() - 2;
		if (unknowns_count == 0) return;

		auto dset = [&](int index, double value)
		{
			values.SetDxy(row_index, index, value);
		};
		auto rget = [&](int index)
		{
			return values.Dx(row_index, index);
		};

		auto dlast = values.Dxy(row_index, values.ColumnsCount() - 1);
		auto dfirst = values.Dxy(row_index, 0);

		SolveTridiagonal(rget, PrecalculatedHY(), dfirst, dlast,
			unknowns_count, dset);
	}

	bool DeBoorKnotsGenerator::IsParallel()
	{
		return is_parallel_;
	}

	const DeBoorKnotsGenerator::Precalculated&
	DeBoorKnotsGenerator::
	PrecalculatedHX() const
	{
		return precalculated_hx_;
	}

	const DeBoorKnotsGenerator::Precalculated&
	DeBoorKnotsGenerator::
	PrecalculatedHY() const
	{
		return precalculated_hy_;
	}

	void 
	DeBoorKnotsGenerator::SetPrecalculatedHX(
		Precalculated precalculated)
	{
		precalculated_hx_ = std::move(precalculated);
	}

	void DeBoorKnotsGenerator::SetPrecalculatedHY(
		Precalculated precalculated)
	{
		precalculated_hy_ = std::move(precalculated);
	}


	void DeBoorKnotsGenerator::InParallel(bool value)
	{
		is_parallel_ = value;
		auto threads = utils::num_threads;
		if (value)
		{
			tridagonals_.reserve(threads);
			for (auto i = tridagonals_.size(); i < threads; i++)
			{
				// create copy of tridiagonal solver
				auto copy_of_first = tridagonals_[0];
				tridagonals_.push_back(std::move(copy_of_first));
			}
		}
		else
		{
			for (int i = 0; i < tridagonals_.size() - 1; ++i) {
				tridagonals_.pop_back();
			}
		}
	}

	void DeBoorKnotsGenerator::InitializeBuffers(const size_t u_count,
		const size_t v_count)
	{
		auto size = std::max(u_count - 2, v_count - 2);
		for (size_t i = 0; i < tridagonals_.size(); i++)
		{
			tridagonals_[i].ResizeBuffers(size);
		}
	}

}