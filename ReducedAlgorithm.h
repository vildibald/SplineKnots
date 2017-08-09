#pragma once

#include "MathFunction.h"
#include "Spline.h"
#include "SurfaceDimension.h"
#include <functional>
#include <vector>
#include <omp.h>

class ReducedAlgorithm {
    splineknots::SurfaceDimension xDimension, yDimension;
    std::vector<double> rightSideX, rightSideY, luBufferX, luBufferY, rightSideXY, luBufferXY,
    rightSideYX, luBufferYX;
    splineknots::InterpolativeMathFunction f;
    int numEquationsX;
    int numEquationsY;

public:
    ReducedAlgorithm(const splineknots::SurfaceDimension xDimension,
                  const splineknots::SurfaceDimension yDimension,
                  const splineknots::InterpolativeMathFunction f);

    splineknots::Spline Calculate(double *time = nullptr);
private:
    void Initialize(splineknots::Spline &spline);

    void FillDx(splineknots::Spline &spline);

    void FillDy(splineknots::Spline &spline);

    void FillDxy(splineknots::Spline &spline);

    void FillDyx(splineknots::Spline &spline);


};


