#include "stdafx.h"
#include "Spline.h"
#include "utils.h"
#include <iostream>

void splineknots::Spline::Print()
{
	using namespace std;
	cout << "---------- Knot matrix ----------" << endl;
	for (size_t i = 0; i < xdim.knotCount; i++)
	{
		cout << "Row " << i << " :\n";
		for (size_t j = 0; j < ydim.knotCount; j++)
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
//
splineknots::Spline::Spline()
	:xdim(0,0,0),ydim(0,0,0),
    z(), dx(), dy(), dxy()
{
}

splineknots::Spline splineknots::Spline::NullMatrix()
{
	Spline nullval;
	return nullval;
}

bool splineknots::Spline::IsNull()
{
    if (xdim.knotCount < 1 || ydim.knotCount < 1)
		return true;
	return false;
}

splineknots::Spline::Spline(SurfaceDimension rowdimension, SurfaceDimension columndimension)
	: xdim(rowdimension), ydim(columndimension)
, z(), dx(), dy(), dxy()
{
	z.resize(xdim.knotCount);
	dx.resize(xdim.knotCount);
	dy.resize(xdim.knotCount);
	dxy.resize(xdim.knotCount);

	for (int i = 0; i < xdim.knotCount; ++i) {
		z[i].resize(ydim.knotCount);
		dx[i].resize(ydim.knotCount);
		dy[i].resize(ydim.knotCount);
		dxy[i].resize(ydim.knotCount);
	}
}