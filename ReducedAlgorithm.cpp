//
// Created by Viliam on 30.4.2017.
//

#include "ReducedAlgorithm.h"
#include "StopWatch.h"
#include "utils.h"

ReducedAlgorithm::ReducedAlgorithm(const splineknots::SurfaceDimension xDimension,
                                   const splineknots::SurfaceDimension yDimension,
                                   const splineknots::InterpolativeMathFunction f)
        : xDimension(xDimension), yDimension(yDimension), rightSideX(), rightSideY(), luBufferX(),
          luBufferY(), rightSideXY(), luBufferXY(), rightSideYX(), luBufferYX(), f(f) {
    numEquationsX = (xDimension.knotCount - 2) / 2;
    numEquationsY = (yDimension.knotCount - 2) / 2;
    rightSideX.resize(numEquationsX);
    rightSideY.resize(numEquationsY);
    luBufferX.resize(numEquationsX);
    luBufferY.resize(numEquationsY);
    rightSideXY.resize(xDimension.knotCount - 2);
    luBufferXY.resize(xDimension.knotCount - 2);
    rightSideYX.resize(yDimension.knotCount - 2);
    luBufferYX.resize(yDimension.knotCount - 2);

}

splineknots::Spline ReducedAlgorithm::Calculate(double *time) {
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

void ReducedAlgorithm::Initialize(splineknots::Spline &spline) {
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

void ReducedAlgorithm::FillDx(splineknots::Spline &spline) {
    double three_div_h = 3 / xDimension.H();
    double twelwe_div_h = 12 / xDimension.H();
    double three_div_4h = 3 / (4 * xDimension.H());
    for (size_t j = 0; j < spline.ColumnsCount(); ++j) {
        for (size_t i = 2, k = 0; i < spline.RowsCount() - 2; i += 2, ++k) {
            rightSideX[k] = three_div_h * (spline.Z(i + 2, j) - spline.Z(i - 2, j))
                            - twelwe_div_h * (spline.Z(i + 1, j) - spline.Z(i - 1, j));
        }
        rightSideX[0] -= spline.Dx(0, j);
        rightSideX[rightSideX.size() - 1] -= spline.Dx(spline.RowsCount() - 1, j);
        utils::SolveDeboorTridiagonalSystemBuffered(-14, &rightSideX.front(), rightSideX.size(),
                                                    &luBufferX.front());
        for (int i = 2, k = 0; i < spline.RowsCount() - 2; i += 2, ++k) {
            spline.SetDx(i, j, rightSideX[k]);
        }
        for (int i = 1; i < spline.RowsCount() - 1; i += 2) {
            spline.SetDx(i, j,
                         three_div_4h * (spline.Z(i + 1, j) - spline.Z(i - 1, j))
                         - 0.25 * (spline.Dx(i + 1, j) + spline.Dx(i - 1, j))
            );
        }
    }
}

void ReducedAlgorithm::FillDy(splineknots::Spline &spline) {
    double three_div_h = 3 / xDimension.H();
    double twelwe_div_h = 12 / xDimension.H();
    double three_div_4h = 3 / (4 * xDimension.H());
    for (size_t i = 0; i < spline.RowsCount(); ++i) {
        for (size_t j = 2, k = 0; j < spline.RowsCount() - 2; j += 2, ++k) {
            rightSideY[k] = three_div_h * (spline.Z(i, j + 2) - spline.Z(i, j - 2))
                            - twelwe_div_h * (spline.Z(i, j + 1) - spline.Z(i, j - 1));
        }
        rightSideY[0] -= spline.Dy(i, 0);
        rightSideY[rightSideY.size() - 1] -= spline.Dy(i, spline.ColumnsCount() - 1);
        utils::SolveDeboorTridiagonalSystemBuffered(-14, &rightSideY.front(), rightSideY.size(),
                                                    &luBufferX.front());
        for (int j = 2, k = 0; j < spline.ColumnsCount() - 2; j += 2, ++k) {
            spline.SetDy(i, j, rightSideY[k]);
        }
        for (int j = 1; j < spline.ColumnsCount() - 1; j += 2) {
            spline.SetDy(i, j,
                         three_div_4h * (spline.Z(i, j + 1) - spline.Z(i, j - 1))
                         - 0.25 * (spline.Dy(i, j + 1) + spline.Dy(i, j - 1))
            );
        }
    }
}

void ReducedAlgorithm::FillDxy(splineknots::Spline &spline) {
    double three_div_h = 3 / xDimension.H();
    int j = 0;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideXY[i - 1] = three_div_h * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideXY[0] -= spline.Dxy(0, j);
    rightSideXY[rightSideXY.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideXY.front(), rightSideXY.size(),
                                                &luBufferXY.front());
    for (int i = 0; i < rightSideXY.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideXY[i]);
    }

    j = spline.ColumnsCount() - 1;
    for (size_t i = 1; i < spline.RowsCount() - 1; ++i) {
        rightSideXY[i - 1] = three_div_h * (spline.Dy(i + 1, j) - spline.Dy(i - 1, j));
    }
    rightSideXY[0] -= spline.Dxy(0, j);
    rightSideXY[rightSideXY.size() - 1] -= spline.Dxy(spline.RowsCount() - 1, j);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideXY.front(), rightSideXY.size(),
                                                &luBufferX.front());
    for (int i = 0; i < rightSideXY.size(); ++i) {
        spline.SetDxy(i + 1, j, rightSideXY[i]);
    }

    ////////////

    three_div_h = 3 / yDimension.H();
    int i = 0;
    for (size_t j = 1; j < spline.ColumnsCount() - j; ++j) {
        rightSideYX[j - 1] = three_div_h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1));
    }
    rightSideYX[0] -= spline.Dxy(i, 0);
    rightSideYX[rightSideYX.size() - 1] -= spline.Dxy(i, spline.ColumnsCount() - 1);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideYX.front(), rightSideYX.size(),
                                                &luBufferYX.front());
    for (int j = 0; j < rightSideYX.size(); ++j) {
        spline.SetDxy(i, j + 1, rightSideYX[j]);
    }

    i = spline.RowsCount() - 1;
    for (size_t j = 1; j < spline.ColumnsCount() - 1; ++j) {
        rightSideYX[j - 1] = three_div_h * (spline.Dy(i, j + 1) - spline.Dy(i, j + 1));
    }
    rightSideYX[0] -= spline.Dxy(i, 0);
    rightSideYX[rightSideYX.size() - 1] -= spline.Dxy(i, spline.ColumnsCount() - 1);
    utils::SolveDeboorTridiagonalSystemBuffered(4, &rightSideYX.front(), rightSideYX.size(),
                                                &luBufferYX.front());
    for (int j = 0; j < rightSideYX.size(); ++j) {
        spline.SetDxy(i, j + 1, rightSideYX[j]);
    }


}

void ReducedAlgorithm::FillDyx(splineknots::Spline &spline) {
    double hx = xDimension.H();
    double hy = yDimension.H();
    double one_div_7 = 1.0 / 7.0;
    double three_div_7hx = 3.0 / (7.0 * hx);
    double three_div_7hy = 3.0 / (7.0 * hy);
    double nine_div_7hx = 9.0 / (7.0 * hx);
    double nine_div_7hxhy = 9.0 / (7.0 * hx * hy);
    double twelve_div_7hx = 12.0 / (7.0 * hx);
    double twelve_div_7hy = 12.0 / (7.0 * hy);
    double three_div_hy = 3.0 / hy;
    double twentyseven_div_7hxhy = 27.0 / (7.0 * hx * hy);
    double thirtysix_div_7hxhy = 36.0 / (7.0 * hx * hy);
    double six_div_hx = 6.0 / hx;
    double twelve_div_hy = 12.0 / hy;
    double onehundredeight_div_7hxhy = 108.0 / (7.0 * hx * hy);
    double eighteen_div_hx = 18.0 / hx;
    double onehundredfortyfour_div_7hxhy = 144.0 / (7.0 * hx * hy);
    double twentyfour_div_hx = 24.0 / hx;
    for (int i = 2; i < spline.RowsCount() - 2; i += 2) {
        int k = 0;
        for (int j = 2; j < spline.ColumnsCount() - 2; j += 2, ++k) {
            rightSideY[k] =
                    one_div_7 * (spline.Dxy(i - 2, j + 2) + spline.Dxy(i - 2, j - 2))
                    - 2 * spline.Dxy(i - 2, j)
                    + three_div_7hx * (spline.Dy(i - 2, j + 2) + spline.Dy(i - 2, j - 2))
                    + three_div_7hy * (-spline.Dx(i - 2, j + 2) + spline.Dx(i - 2, j - 2))
                    + nine_div_7hx * (spline.Dy(i, j + 2) + spline.Dy(i, j - 2))
                    + nine_div_7hxhy * (-spline.Z(i - 2, j + 2) + spline.Z(i - 2, j - 2))
                    + twelve_div_7hx * (-spline.Dy(i - 1, j + 2) - spline.Dy(i - 1, j - 2))
                    + twelve_div_7hy * (spline.Dx(i - 2, j + 1) - spline.Dx(i - 2, j - 1))
                    + three_div_hy * (spline.Dx(i, j + 2) - spline.Dx(i, j - 2))
                    + twentyseven_div_7hxhy * (-spline.Z(i, j + 2) + spline.Z(i, j - 2))
                    + thirtysix_div_7hxhy * (spline.Z(i - 1, j + 2) - spline.Z(i - 1, j - 2)
                                             + spline.Z(i - 2, j + 1) - spline.Z(i - 2, j - 1))
                    - six_div_hx * spline.Dy(i - 2, j)
                    + twelve_div_hy * (spline.Dx(i, j + 1) + spline.Dx(i, j - 1))
                    + onehundredeight_div_7hxhy * (spline.Z(i, j + 1) - spline.Z(i, j - 1))
                    - eighteen_div_hx * spline.Dy(i, j)
                    + onehundredfortyfour_div_7hxhy * (-spline.Z(i - 1, j + 1)
                                                       + spline.Z(i - 1, j - 1))
                    + twentyfour_div_hx * spline.Dy(i - 1, j);

        }
        rightSideY[0] -= spline.Dxy(i, 0);
        rightSideY[k - 1] -= spline.Dxy(i, spline.ColumnsCount() - 1);
        utils::SolveDeboorTridiagonalSystemBuffered(-14, &rightSideY.front(), k,
                                                    &luBufferY.front());
        for (int j = 2, k = 0; j < spline.ColumnsCount() - 2; j += 2, ++k) {
            spline.SetDxy(i, j, rightSideY[k]);
        }
    }

    double one_div_16 = 1.0 / 16.0;
    double three_div_16hy = 3.0 / (16.0 * hy);
    double three_div_16hx = 3.0 / (16.0 * hx);
    double nine_div_16hxhy = 9.0 / (16.0 * hx * hy);
    for (int i = 1; i < spline.RowsCount() - 1; i += 2) {
        for (int j = 1; j < spline.ColumnsCount() - 1; j += 2) {
            spline.SetDxy(i, j,
                          one_div_16 * (spline.Dxy(i + 1, j + 1) + spline.Dxy(i + 1, j - 1)
                                        + spline.Dxy(i - 1, j + 1) + spline.Dxy(i - 1, j - 1))
                          - three_div_16hy * (spline.Dx(i + 1, j + 1) - spline.Dx(i + 1, j - 1)
                                              + spline.Dx(i - 1, j + 1) - spline.Dx(i - 1, j - 1))
                          - three_div_16hx * (spline.Dy(i + 1, j + 1) + spline.Dy(i + 1, j - 1)
                                              - spline.Dy(i - 1, j + 1) - spline.Dy(i - 1, j - 1))
                          + nine_div_16hxhy * (spline.Z(i + 1, j + 1) - spline.Z(i + 1, j - 1)
                                               - spline.Z(i - 1, j + 1) + spline.Z(i - 1, j - 1))
            );
        }
    }

    double three_div_4h = 3 / (4 * xDimension.H());
//    for (int j = 1; j < spline.RowsCount() - 1; j += 2) {
//        for (int i = 2; i < spline.ColumnsCount() - 2; i+=2) {
//            spline.SetDxy(i, j,
//                         three_div_4h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1))
//                         - 0.25 * (spline.Dxy(i, j + 1) + spline.Dxy(i,j - 1))
//            );
//        }
//        for (int i = 1; i < spline.ColumnsCount() - 1; i+=2) {
//            spline.SetDxy(i, i,
//                          three_div_4h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1))
//                          - 0.25 * (spline.Dxy(i, j + 1) + spline.Dxy(i, j - 1))
//            );
//        }
//    }
    for (int i = 1; i < spline.RowsCount() - 1; i+=2) {
        for (int j = 2; j < spline.ColumnsCount() - 2; j+=2) {
            spline.SetDxy(i, j,
                         three_div_4h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1))
                         - 0.25 * (spline.Dxy(i, j + 1) + spline.Dxy(i,j - 1))
            );
        }
    }
    for (int i = 2; i < spline.RowsCount() - 2; i+=2) {
        for (int j = 1; j < spline.ColumnsCount() - 1; j+=2) {
            spline.SetDxy(i, j,
                          three_div_4h * (spline.Dx(i, j + 1) - spline.Dx(i, j - 1))
                          - 0.25 * (spline.Dxy(i, j + 1) + spline.Dxy(i,j - 1))
            );
        }
    }
}

