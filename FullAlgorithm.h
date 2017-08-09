#pragma once

#include "MathFunction.h"
#include "Spline.h"
#include "SurfaceDimension.h"
#include <functional>
#include <vector>
#include <omp.h>

class FullAlgorithm {
    splineknots::SurfaceDimension xDimension, yDimension;
    std::vector<double> rightSideX, rightSideY, luBufferX, luBufferY;
    splineknots::InterpolativeMathFunction f;
    int numEquationsX;
    int numEquationsY;

public:
    FullAlgorithm(const splineknots::SurfaceDimension xDimension,
                  const splineknots::SurfaceDimension yDimension,
                  const splineknots::InterpolativeMathFunction f);

    splineknots::Spline Calculate(double *time = nullptr);

private:
    void Initialize(splineknots::Spline &values);

    void FillDx(splineknots::Spline &values);

    void FillDy(splineknots::Spline &values);

    void FillDxy(splineknots::Spline &values);

    void FillDyx(splineknots::Spline &values);
};
