#include "stdafx.h"
#include "SplineKnots.h"
#include <iostream>
#include <algorithm>
#include "OperationsBenchmark.h"
#include "ComparisonBenchmarkResult.h"
#include <numeric>
#include "StopWatch.h"
#include "FullAlgorithm.h"
#include "ReducedAlgorithm.h"

//
//void EqualityComparison();
//
void MulDivBenchmark() {
    OperationsBenchmark bencher;
    bencher.BenchAll();
}


void PrintSurfaceDeboorResult(ComparisonBenchmarkResult result) {
    std::cout << "Full : " << result.FirstAlg() << std::endl;
    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
    std::cout << "Difference F/R: " << result.Ratio() << std::endl;
}

void SurfaceBenchmark(int numIterations, int numKnots) {
    splineknots::SurfaceDimension uDim(-20, 20, numKnots);
    splineknots::SurfaceDimension vDim = uDim;

    splineknots::MathFunction function = [](double x, double y) {
        return sin(sqrt(x * x + y * y));
    };
    splineknots::InterpolativeMathFunction f = function;

    FullAlgorithm full(uDim, vDim, f);
    ReducedAlgorithm reduced(uDim, vDim, f);

    std::vector<double> calculatedResults;
    std::vector<double> fullTimes;
    fullTimes.reserve(numIterations);
    std::vector<double> reducedTimes;
    reducedTimes.reserve(numIterations);

    calculatedResults.reserve(numIterations * 2);

    for (size_t i = 0; i < numIterations; i++) {
        double time = 0;
        auto result = full.Calculate(&time);
        calculatedResults.emplace_back(result.Dxy(1, 1));
        fullTimes.emplace_back(time);
    }

    for (size_t i = 0; i < numIterations; i++) {
        double time = 0;
        auto result = reduced.Calculate(&time);
        calculatedResults.emplace_back(result.Dxy(1, 1));
        reducedTimes.emplace_back(time);
    }

    auto full_time = static_cast<double>(std::accumulate(fullTimes.begin(),
                                                         fullTimes.end(), 0))
                     / static_cast<double>(numIterations);
    auto reduced_time = static_cast<double>(std::accumulate(
            reducedTimes.begin(), reducedTimes.end(), 0))
                        / static_cast<double>(numIterations);
    std::cout << "Ignore " << calculatedResults[0] << std::endl;
    PrintSurfaceDeboorResult(ComparisonBenchmarkResult(full_time, reduced_time));
}

void EqualityComparison() {
    splineknots::SurfaceDimension udimension(-3, 3, 7);
    splineknots::SurfaceDimension vdimension = udimension;

    splineknots::MathFunction function = [](double x, double y) {
        return sin(sqrt(x * x + y * y));
    };
    splineknots::InterpolativeMathFunction f = function;

    FullAlgorithm full(udimension, vdimension, f);
    ReducedAlgorithm reduced(udimension, vdimension, f);

    std::vector<double> calculatedResults;
    std::vector<double> fullTimes;
    std::vector<double> reducedTimes;

    auto resultFull = full.Calculate();
    auto resultReduced = reduced.Calculate();

    splineknots::InterpolativeMathFunction interpolativeMathFunction = function;
    std::cout << "---------- Knot matrix ----------" << std::endl;

    auto minDiffDx = std::numeric_limits<double>::max(), maxDiffDx = std::numeric_limits<double>::min();
    auto minDiffDy = std::numeric_limits<double>::max(), maxDiffDy = std::numeric_limits<double>::min();
    auto minDiffDxy = std::numeric_limits<double>::max(), maxDiffDxy = std::numeric_limits<double>::min();

    const double hx = udimension.H();
    const double hy = vdimension.H();

    double x = udimension.min;
    for (unsigned int i = 0; i < resultFull.RowsCount(); ++i, x += hx) {
//    for (unsigned int i = numKnots/2; i < numKnots/2+1; ++i) {
//    for (unsigned int i = 4; i < 5; ++i) {
        std::cout << "Row " << i << " :\n";
        double y = vdimension.min;
        for (unsigned int j = 0; j < resultFull.ColumnsCount(); ++j, y += hy) {
            minDiffDx = std::min(minDiffDx,
                                 std::abs(resultFull.Dx(i, j) - resultReduced.Dx(i, j)));
            minDiffDy = std::min(minDiffDy,
                                 std::abs(resultFull.Dy(i, j) - resultReduced.Dy(i, j)));
            minDiffDxy = std::min(minDiffDxy,
                                  std::abs(resultFull.Dxy(i, j) - resultReduced.Dxy(i, j)));

            maxDiffDx = std::max(maxDiffDx,
                                 std::abs(resultFull.Dx(i, j) - resultReduced.Dx(i, j)));
            maxDiffDy = std::max(maxDiffDy,
                                 std::abs(resultFull.Dy(i, j) - resultReduced.Dy(i, j)));
            maxDiffDxy = std::max(maxDiffDxy,
                                  std::abs(resultFull.Dxy(i, j) - resultReduced.Dxy(i, j)));

            std::cout << "\nColumn " << j << ":\n"
                      << "M. z: " << interpolativeMathFunction.Z()(x, y) << '\n'
                      << "M. dx: " << interpolativeMathFunction.Dx()(x, y) << '\n'
                      << "M. dy: " << interpolativeMathFunction.Dy()(x, y) << '\n'
                      << "M. dxy: " << interpolativeMathFunction.Dxy()(x, y) << "\n\n";
            std::cout << "F. z: " << resultFull.Z(i, j) << '\n'
                      << "F. dx: " << resultFull.Dx(i, j) << '\n'
                      << "F. dy: " << resultFull.Dy(i, j) << '\n'
                      << "F. dxy: " << resultFull.Dxy(i, j) << "\n\n";
            std::cout << "R. z: " << resultReduced.Z(i, j) << '\n'
                      << "R. dx: " << resultReduced.Dx(i, j) << '\n'
                      << "R. dy: " << resultReduced.Dy(i, j) << '\n'
                      << "R. dxy: " << resultReduced.Dxy(i, j) << '\n';
        }
        std::cout << std::endl;
    }

    std::cout << "-------------------------------" << std::endl;

    std::cout << "Min diff Dx: " << minDiffDx << std::endl;
    std::cout << "Max diff Dx: " << maxDiffDx << std::endl;
    std::cout << "Min diff Dy: " << minDiffDy << std::endl;
    std::cout << "Max diff Dy: " << maxDiffDy << std::endl;
    std::cout << "Min diff Dxy: " << minDiffDxy << std::endl;
    std::cout << "Max diff Dxy: " << maxDiffDxy << std::endl;
}


using namespace splineknots;


int main() {

    while (true) {
        //std::cout << clock();
        // Console clear ...
        // ... for Windows,
        system("cls");
        // ... for Linux/Unix.
        //system("clear");
        std::cout << "1: Instructions benchmark." << std::endl;
        std::cout << "2: Spline surface benchmark." << std::endl;
        std::cout << "3: Equality comparison." << std::endl;
        std::cout << "Q: End program" << std::endl;
        char input;
        std::cin >> input;
        std::cin.get();
        std::cout << std::endl << "---------------" << std::endl;
        unsigned int num_iterations;
        unsigned int num_knots;
        switch (input) {
            case '1':
                std::cout << "Instructions benchmark" << std::endl << std::endl;
                MulDivBenchmark();
                break;
            case '2':
                std::cout << "Spline surface benchmark" << std::endl << std::endl;
                num_iterations = 10;
                num_knots = 2001;
                SurfaceBenchmark(num_iterations, num_knots);
                break;
            case '3':
                EqualityComparison();
                break;
            case 'q':
            case 'Q':
                return 0;
        }

        std::cout << "===================" << std::endl;
        system("pause");
    }
    return 0;

}

