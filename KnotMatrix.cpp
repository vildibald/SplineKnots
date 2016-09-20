#include "stdafx.h"
#include "KnotMatrix.h"
#include "utils.h"
#include <iostream>

//#include "utils_template.cpp"


splineknots::KnotMatrix::~KnotMatrix() noexcept
{
	if (!z_) return;
	for (size_t i = 0; i < xdim_.knot_count; i++)
	{
		delete[] z_[i];
		delete[] dy_[i];
		delete[] dxy_[i];
	}
	for (size_t i = 0; i < ydim_.knot_count; i++)
	{
		delete[] dx_[i];
	}
	delete[] z_;
	delete[] dx_;
	delete[] dy_;
	delete[] dxy_;
	z_ = nullptr;
	dx_ = nullptr;
	dy_ = nullptr;
	dxy_ = nullptr;
}

void splineknots::KnotMatrix::Print()
{
	using namespace std;
	cout << "---------- Knot matrix ----------" << endl;
	for (size_t i = 0; i < xdim_.knot_count; i++)
	{
		cout << "Row " << i << " :\n";
		for (size_t j = 0; j < ydim_.knot_count; j++)
		{
			cout << j << ":\n"
				<< "z: " << z_[i][j] << '\n'
				<< "dx: " << dx_[j][i] << '\n'
				<< "dy: " << dy_[i][j] << '\n'
				<< "dxy: " << dxy_[i][j] << '\n';
		}
		cout << endl;
	}
	cout << "-------------------------------" << endl;
}

splineknots::KnotMatrix::KnotMatrix()
	:xdim_(0,0,0),ydim_(0,0,0),
	 z_(nullptr),
	 dx_(nullptr), dy_(nullptr), dxy_(nullptr)
{
}

splineknots::KnotMatrix splineknots::KnotMatrix::NullMatrix()
{
	KnotMatrix nullval;
	return nullval;
}

bool splineknots::KnotMatrix::IsNull()
{
	if (!z_ || !dx_ || !dy_ || !dxy_ || xdim_.knot_count < 1 || ydim_.knot_count < 1)
		return true;
	return false;
}

splineknots::KnotMatrix::KnotMatrix(SurfaceDimension rowdimension, SurfaceDimension columndimension)
	: xdim_(std::move(rowdimension)), ydim_(std::move(columndimension))
{
	z_ = new double*[xdim_.knot_count];
	dy_ = new double*[xdim_.knot_count];
	dxy_ = new double*[xdim_.knot_count];
	for (size_t i = 0; i < xdim_.knot_count; i++)
	{
		z_[i] = new double[ydim_.knot_count];
		dy_[i] = new double[ydim_.knot_count];
		dxy_[i] = new double[ydim_.knot_count];
	}
	dx_ = new double*[ydim_.knot_count];
	for (size_t i = 0; i < ydim_.knot_count; i++)
	{
		dx_[i] = new double[xdim_.knot_count];
	}
}

splineknots::KnotMatrix::KnotMatrix(const KnotMatrix& other)
	: xdim_(other.xdim_),
	ydim_(other.ydim_)
{
	z_ = new double*[xdim_.knot_count];
	dy_ = new double*[xdim_.knot_count];
	dxy_ = new double*[xdim_.knot_count];
	for (size_t i = 0; i < other.RowsCount(); i++)
	{
		z_[i] = new double[ydim_.knot_count];
		dy_[i] = new double[ydim_.knot_count];
		dxy_[i] = new double[ydim_.knot_count];
		memcpy(z_[i], other.z_[i], ydim_.knot_count);
		memcpy(dy_[i], other.dy_[i], ydim_.knot_count);
		memcpy(dxy_[i], other.dxy_[i], ydim_.knot_count);
	}
	dx_ = new double*[ydim_.knot_count];
	for (size_t i = 0; i < other.ColumnsCount(); i++)
	{
		dx_[i] = new double[xdim_.knot_count];
		memcpy(dx_[i], other.dx_[i], xdim_.knot_count);
	}
}

splineknots::KnotMatrix::KnotMatrix(KnotMatrix&& other) noexcept
	: xdim_(std::move(other.xdim_)),
	ydim_(std::move(other.ydim_))
{
	z_ = other.z_;
	other.z_ = nullptr;
	dx_ = other.dx_;
	other.dx_ = nullptr;
	dy_ = other.dy_;
	other.dy_ = nullptr;
	dxy_ = other.dxy_;
	other.dxy_ = nullptr;
}

splineknots::KnotMatrix& splineknots::KnotMatrix::operator=(const KnotMatrix&
	other)
{
	if (&other != this)
	{
		utils::DeleteJaggedArray(z_, xdim_.knot_count, ydim_.knot_count);
		utils::DeleteJaggedArray(dx_, ydim_.knot_count, xdim_.knot_count);
		utils::DeleteJaggedArray(dy_, xdim_.knot_count, ydim_.knot_count);
		utils::DeleteJaggedArray(dxy_, xdim_.knot_count, ydim_.knot_count);
		xdim_ = other.xdim_;
		ydim_ = other.ydim_;
		z_ = utils::CreateJaggedArray<double>(xdim_.knot_count, ydim_.knot_count);
		dx_ = utils::CreateJaggedArray<double>(ydim_.knot_count, xdim_.knot_count);
		dy_ = utils::CreateJaggedArray<double>(xdim_.knot_count, ydim_.knot_count);
		dxy_ = utils::CreateJaggedArray<double>(xdim_.knot_count, ydim_.knot_count);
		for (size_t i = 0; i < xdim_.knot_count; i++)
		{
			memcpy(z_[i], other.z_[i], ydim_.knot_count);
			memcpy(dy_[i], other.dy_[i], ydim_.knot_count);
			memcpy(dxy_[i], other.dxy_[i], ydim_.knot_count);
		}
		for (size_t i = 0; i < ydim_.knot_count; i++)
		{
			memcpy(dx_[i], other.dx_[i], xdim_.knot_count);
		}
	}
	return *this;
}

splineknots::KnotMatrix& splineknots::KnotMatrix::operator=(KnotMatrix&& other) noexcept
{
	if (&other != this)
	{
		utils::DeleteJaggedArray(z_, xdim_.knot_count, ydim_.knot_count);
		utils::DeleteJaggedArray(dx_, ydim_.knot_count, xdim_.knot_count);
		utils::DeleteJaggedArray(dy_, xdim_.knot_count, ydim_.knot_count);
		utils::DeleteJaggedArray(dxy_, xdim_.knot_count, ydim_.knot_count);
		xdim_ = std::move(other.xdim_);
		ydim_ = std::move(other.ydim_);
		z_ = other.z_;
		other.z_ = nullptr;
		dx_ = other.dx_;
		other.dx_ = nullptr;
		dy_ = other.dy_;
		other.dy_ = nullptr;
		dxy_ = other.dxy_;
		other.dxy_ = nullptr;

	}
	return *this;
}
