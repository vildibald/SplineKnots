#include "FullAlgorithm.h"
#include "StopWatch.h"
#include "utils.h"

FullAlgorithm::FullAlgorithm(const splineknots::SurfaceDimension xDimension,
                             const splineknots::SurfaceDimension yDimension,
                             const splineknots::InterpolativeMathFunction f)
        : xDimension(xDimension), yDimension(yDimension), rightSideX(), rightSideY(), luBufferX(),
          luBufferY(), f(f) {
    numEquationsX = xDimension.knotCount - 2;
    numEquationsY = yDimension.knotCount - 2;
    rightSideX.resize(numEquationsX);
    rightSideY.resize(numEquationsY);
    luBufferX.resize(numEquationsX);
    luBufferY.resize(numEquationsY);

}

splineknots::Spline FullAlgorithm::Calculate(double *time) {
    StopWatch sw;
    using namespace splineknots;
    Spline spline{xDimension, yDimension};
    Initialize(spline);

    sw.Start();
    FillDx(spline);
    FillDy(spline);
    FillDxy(spline);
    FillDyx(spline);
    sw.Stop();

    if (time != nullptr) {
        *time = sw.EllapsedTime();
    }
    return spline;
}


void FullAlgorithm::FillDx(splineknots::Spline &spline) {
    double three_div_h = 3 / xDimension.H();
    for (size_t j = 0; j < spline.ColumnsCount(); ++j) {
        for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
            rightSideX[i - 1] = three_div_h * (spline.Z(i + 1, j) - spline.Z(i - 1, j));
        }
        rightSideX[0] -= spline.Dx(0, j);
        rightSideX[rightSideX.size() - 1] -= spline.Dx(spline.RowsCount() - 1, j);
        utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideX.front(), rightSideX.size(),
                                                    &luBufferX.front());
        for (int i = 0; i < rightSideX.size(); ++i) {
            spline.SetDx(i + 1, j, rightSideX[i]);
        }
    }
}

void FullAlgorithm::FillDy(splineknots::Spline &spline) {
    double three_div_h = 3 / yDimension.H();
    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        for (size_t j = 1; j < spline.RowsCount() - 1; ++j) {
            rightSideY[j - 1] = three_div_h * (spline.Z(i, j + 1) - spline.Z(i, j - 1));
        }
        rightSideY[0] -= spline.Dy(i, 0);
        rightSideY[rightSideY.size() - 1] -= spline.Dy(i, spline.ColumnsCount() - 1);
        utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideY.front(), rightSideY.size(),
                                                    &luBufferY.front());
        for (int j = 0; j < rightSideY.size(); ++j) {
            spline.SetDy(i, j + 1, rightSideY[j]);
        }
    }
}

void FullAlgorithm::FillDxy(splineknots::Spline &spline) {
    double three_div_h = 3 / xDimension.H();
    int j = 0;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideX[i - 1] = three_div_h * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideX[0] -= spline.Dxy(0, j);
    rightSideX[rightSideX.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideX.front(), rightSideX.size(),
                                                &luBufferX.front());
    for (int i = 0; i < rightSideX.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideX[i]);
    }

    j = spline.ColumnsCount() - 1;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideX[i - 1] = three_div_h * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideX[0] -= spline.Dxy(0, j);
    rightSideX[rightSideX.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideX.front(), rightSideX.size(),
                                                &luBufferX.front());
    for (int i = 0; i < rightSideX.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideX[i]);
    }
}

void FullAlgorithm::FillDyx(splineknots::Spline &spline) {
    double three_div_h = 3 / yDimension.H();
    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        for (size_t j = 1; j < spline.RowsCount() - 1; ++j) {
            rightSideY[j - 1] = three_div_h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1));
        }
        rightSideY[0] -= spline.Dxy(i, 0);
        rightSideY[rightSideY.size() - 1] -= spline.Dxy(i, spline.ColumnsCount() - 1);
        utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideY.front(), rightSideY.size(),
                                                    &luBufferX.front());
        for (int j = 0; j < rightSideY.size(); ++j) {
            spline.SetDxy(i, j + 1, rightSideY[j]);
        }
    }
}

void FullAlgorithm::Initialize(splineknots::Spline &spline) {
    auto hx = xDimension.H();
    auto hy = yDimension.H();
    auto u = xDimension.min;
    // Init Z
    for (auto i = 0; i < xDimension.knotCount; i++, u += hx) {
        auto v = xDimension.min;
        for (auto j = 0; j < yDimension.knotCount; j++, v += hy) {
            auto z = f.Z()(u, v);
            spline.SetZ(i, j, z);
        }
    }
    // Init Dx
    u = xDimension.min;
    auto uKnotCountMin1 = xDimension.knotCount - 1;
    for (auto j = 0; j < yDimension.knotCount; j++) {
//            auto v = vDimension.min;
        auto x0 = spline.X(0);
        auto y0 = spline.Y(0);
        auto x1 = spline.X(uKnotCountMin1);
        auto y1 = spline.Y(uKnotCountMin1);
        auto d0 = f.Dx()(x0, y0);
        auto d1 = f.Dx()(x1, y1);
        spline.SetDx(0, j, d0);
        spline.SetDx(uKnotCountMin1, j, d1);
    }
    // Init Dy
    auto vKnotCountMin1 = yDimension.knotCount - 1;
    for (auto i = 0; i < xDimension.knotCount; i++) {
        auto x0 = spline.X(0);
        auto y0 = spline.Y(0);
        auto x1 = spline.X(vKnotCountMin1);
        auto y1 = spline.Y(vKnotCountMin1);
        auto d0 = f.Dy()(x0, y0);
        auto d1 = f.Dy()(x1, y1);
        spline.SetDy(i, 0, d0);
        spline.SetDy(i, vKnotCountMin1, d1);

    }
    // Init Dxy
    spline.SetDxy(0, 0, f.Dxy()(spline.X(0), spline.Y(0)));
    spline.SetDxy(uKnotCountMin1, 0, f.Dxy()(spline.X(0),
                                             spline.Y(uKnotCountMin1)));
    spline.SetDxy(0, vKnotCountMin1, f.Dxy()(spline.X(vKnotCountMin1),
                                             spline.Y(0)));
    spline.SetDxy(uKnotCountMin1, vKnotCountMin1, f.Dxy()(spline.X(vKnotCountMin1),
                                                          spline.Y(uKnotCountMin1)));
}
    