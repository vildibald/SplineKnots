#include "stdafx.h"
#include "KnotMatrix.h"
#include "utils.h"
#include <iostream>

//#include "utils_template.cpp"


splineknots::KnotMatrix::~KnotMatrix() noexcept
{
	if (!z) return;
	for (size_t i = 0; i < x.size(); i++)
	{
		delete[] z[i];
		delete[] dy[i];
		delete[] dxy[i];
	}
	for (size_t i = 0; i < y.size(); i++)
	{
		delete[] dx[i];
	}
	delete[] z;
	delete[] dx;
	delete[] dy;
	delete[] dxy;
	z = nullptr;
	dx = nullptr;
	dy = nullptr;
	dxy = nullptr;
}

void splineknots::KnotMatrix::Print()
{
	using namespace std;
	cout << "---------- Knot matrix ----------" << endl;
	for (size_t i = 0; i < RowsCount(); i++)
	{
		cout << "Row " << i << " :\n";
		for (size_t j = 0; j < ColumnsCount(); j++)
		{
			cout << j << ":\n"
				<< "z: " << z[i][j] << '\n'
				<< "dx: " << dx[j][i] << '\n'
				<< "dy: " << dy[i][j] << '\n'
				<< "dxy: " << dxy[i][j] << '\n';
		}
		cout << endl;
	}
	cout << "-------------------------------" << endl;
}

splineknots::KnotMatrix::KnotMatrix()
		: x(), y(),
	 z(nullptr),
	 dx(nullptr), dy(nullptr), dxy(nullptr)
{
}

splineknots::KnotMatrix splineknots::KnotMatrix::NullMatrix()
{
	KnotMatrix nullval;
	return nullval;
}

bool splineknots::KnotMatrix::IsNull()
{
	if (!z || !dx || !dy || !dxy || RowsCount() < 1 || ColumnsCount() < 1)
		return true;
	return false;
}

splineknots::KnotMatrix::KnotMatrix(KnotVector rowVector, KnotVector columnVector)
	: x(std::move(rowVector)), y(std::move(columnVector))
{
	z = new double*[RowsCount()];
	dy = new double*[RowsCount()];
	dxy = new double*[RowsCount()];
	for (size_t i = 0; i < RowsCount(); i++)
	{
		z[i] = new double[ColumnsCount()];
		dy[i] = new double[ColumnsCount()];
		dxy[i] = new double[ColumnsCount()];
	}
	dx = new double*[ColumnsCount()];
	for (size_t i = 0; i < ColumnsCount(); i++)
	{
		dx[i] = new double[RowsCount()];
	}
}

splineknots::KnotMatrix::KnotMatrix(const KnotMatrix& other)
	: x(other.x),
	y(other.y)
{
	z = new double*[RowsCount()];
	dy = new double*[RowsCount()];
	dxy = new double*[RowsCount()];
	for (size_t i = 0; i < other.RowsCount(); i++)
	{
		z[i] = new double[ColumnsCount()];
		dy[i] = new double[ColumnsCount()];
		dxy[i] = new double[ColumnsCount()];
		memcpy(z[i], other.z[i], ColumnsCount());
		memcpy(dy[i], other.dy[i], ColumnsCount());
		memcpy(dxy[i], other.dxy[i], ColumnsCount());
	}
	dx = new double*[ColumnsCount()];
	for (size_t i = 0; i < other.ColumnsCount(); i++)
	{
		dx[i] = new double[RowsCount()];
		memcpy(dx[i], other.dx[i], RowsCount());
	}
}

splineknots::KnotMatrix::KnotMatrix(KnotMatrix&& other) noexcept
	: x(std::move(other.x)),
	y(std::move(other.y))
{
	z = other.z;
	other.z = nullptr;
	dx = other.dx;
	other.dx = nullptr;
	dy = other.dy;
	other.dy = nullptr;
	dxy = other.dxy;
	other.dxy = nullptr;
}

splineknots::KnotMatrix& splineknots::KnotMatrix::operator=(const KnotMatrix&
	other)
{
	if (&other != this)
	{
		utils::DeleteJaggedArray(z, RowsCount(), ColumnsCount());
		utils::DeleteJaggedArray(dx, ColumnsCount(), RowsCount());
		utils::DeleteJaggedArray(dy, RowsCount(), ColumnsCount());
		utils::DeleteJaggedArray(dxy, RowsCount(), ColumnsCount());
		x = other.x;
		y = other.y;
		z = utils::CreateJaggedArray<double>(RowsCount(), ColumnsCount());
		dx = utils::CreateJaggedArray<double>(ColumnsCount(), RowsCount());
		dy = utils::CreateJaggedArray<double>(RowsCount(), ColumnsCount());
		dxy = utils::CreateJaggedArray<double>(RowsCount(), ColumnsCount());
		for (size_t i = 0; i < RowsCount(); i++)
		{
			memcpy(z[i], other.z[i], ColumnsCount());
			memcpy(dy[i], other.dy[i], ColumnsCount());
			memcpy(dxy[i], other.dxy[i], ColumnsCount());
		}
		for (size_t i = 0; i < ColumnsCount(); i++)
		{
			memcpy(dx[i], other.dx[i], RowsCount());
		}
	}
	return *this;
}

splineknots::KnotMatrix& splineknots::KnotMatrix::operator=(KnotMatrix&& other) noexcept
{
	if (&other != this)
	{
		utils::DeleteJaggedArray(z, RowsCount(), ColumnsCount());
		utils::DeleteJaggedArray(dx, ColumnsCount(), RowsCount());
		utils::DeleteJaggedArray(dy, RowsCount(), ColumnsCount());
		utils::DeleteJaggedArray(dxy, RowsCount(), ColumnsCount());
		x = std::move(other.x);
		y = std::move(other.y);
		z = other.z;
		other.z = nullptr;
		dx = other.dx;
		other.dx = nullptr;
		dy = other.dy;
		other.dy = nullptr;
		dxy = other.dxy;
		other.dxy = nullptr;

	}
	return *this;
}
