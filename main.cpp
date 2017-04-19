#include "stdafx.h"
#include "SplineKnots.h"
#include <iostream>
#include <algorithm>
#include "MulVsDiv.h"
#include "ComparisonBenchmarkResult.h"
#include <numeric>
#include "StopWatch.h"
#include "CurveDeboorKnotsGenerator.h"
#include "ReducedCurveDeboorGenerator.h"
#include "EnhancedKnotsGenerator.h"

void LUComparison()
{
    using namespace std;
    using namespace utils;
    const size_t num_equations = 9;
    vector<double> right_side_cs(num_equations);
    for (size_t i = 0; i < num_equations; i++)
    {
        right_side_cs[i] = sin(i);
    }
    /*auto dfirst = cos(0);
    auto dlast = cos(num_equations - 1);*/
    auto right_side_vk = right_side_cs;
    auto result_cs = SolveCsabaDeboorTridiagonalSystem(4,
                                                       &right_side_cs.front(), num_equations);
    vector<double> buffer(num_equations);
    SolveDeboorTridiagonalSystemBuffered(1, 4, 1, &right_side_vk.front(),
                                         num_equations, &buffer.front());
    auto result_vk(move(right_side_vk));
    for (size_t i = 0; i < num_equations; i++)
    {
        cout << "Csaba LU: " << result_cs[i] << "; Vilo LU: " << result_vk[i]
             << endl;
    }
    cout << endl;
}

void MulDivBenchmark()
{
    MulVsDiv bencher;
    bencher.BenchAll();
}


ComparisonBenchmarkResult CurveBenchmark(int num_iterations, int num_knots,
                                         bool optimized_lu = true)
{
    const int num_repetitions = 100;
    splineknots::MathFunction function = [](double x, double y)
    {
        return sin(sqrt(x * x));
    };

    splineknots::CurveDeboorKnotsGenerator full(function, optimized_lu);
    splineknots::ReducedCurveDeboorKnotsGenerator reduced(function, optimized_lu);

    splineknots::SurfaceDimension udimension(-2, 2, num_knots);

    std::vector<double> calculated_results;

    std::vector<unsigned int> full_times;
    full_times.reserve(num_iterations);
    std::vector<unsigned int> reduced_times;
    reduced_times.reserve(num_iterations);
    calculated_results.reserve(num_iterations * 2 * num_repetitions);

    for (size_t i = 0; i < num_iterations; i++)
    {
        double time = 0;
        auto result = reduced.GenerateKnots(udimension, &time);
        calculated_results.push_back(result[0]);
        reduced_times.push_back(time);
    }

    for (size_t i = 0; i < num_iterations; i++)
    {
        double time = 0;
        auto result = full.GenerateKnots(udimension, &time);
        calculated_results.push_back(result[0]);
        full_times.push_back(time);
    }

    auto full_time = static_cast<double>(std::accumulate(full_times.begin(), full_times.end(), 0))
                     /num_iterations;
    auto reduced_time = static_cast<double>(std::accumulate(reduced_times.begin(),
                                                            reduced_times.end(), 0))
                        /num_iterations;
    std::cout << "Ignore " << calculated_results[0] << std::endl;
    return ComparisonBenchmarkResult(full_time, reduced_time);
}

ComparisonBenchmarkResult SurfaceBenchmark(int num_iterations, int num_knots,
                                           bool in_parallel = false, bool optimized_lu = true)
{
    splineknots::MathFunction function = [](double x, double y)
    {
        return sin(sqrt(x * x + y * y));
    };

    splineknots::FullKnotsGenerator full(function, optimized_lu);
    splineknots::ReducedKnotsGenerator reduced(function, optimized_lu);
    splineknots::EnhancedKnotsGenerator enhanced(function,
                                                         optimized_lu);
    full.InParallel(in_parallel);
    reduced.InParallel(in_parallel);
    enhanced.InParallel(in_parallel);
    num_knots = num_knots % 2 == 0 ? num_knots + 1 : num_knots;
    const splineknots::SurfaceDimension udimension(-3, 3, num_knots);
    const splineknots::SurfaceDimension vdimension(udimension);

    std::vector<double> calculated_results;
    std::vector<double> full_times;
    full_times.reserve(num_iterations);
    std::vector<double> reduced_times;
    reduced_times.reserve(num_iterations);
    std::vector<double> enhanced_times;
    enhanced_times.reserve(num_iterations);
    calculated_results.reserve(num_iterations * 3);

    for (size_t i = 0; i < num_iterations; i++)
    {
        double time = 0;
        auto result = full.GenerateKnots(udimension, vdimension, &time);
        result.Print();
        calculated_results.push_back(result.Dxy(1, 1));
        full_times.push_back(time);
    }

    for (size_t i = 0; i < num_iterations; i++)
    {
        double time = 0;
        auto result = reduced.GenerateKnots(udimension, vdimension,&time);
        result.Print();
        calculated_results.push_back(result.Dxy(1, 1));
        reduced_times.push_back(time);
    }

    for (size_t i = 0; i < num_iterations; i++)
    {
        double time = 0;
        auto result = enhanced.GenerateKnots(udimension, vdimension,&time);
        result.Print();
        calculated_results.push_back(result.Dxy(1, 1));
        enhanced_times.push_back(time);
    }

    auto full_time = static_cast<double>(std::accumulate(full_times.begin(),
                                                         full_times.end(), 0))
                     / static_cast<double>(num_iterations);
    auto reduced_time = static_cast<double>(std::accumulate(
            reduced_times.begin(), reduced_times.end(), 0))
                        / static_cast<double>(num_iterations);
    auto enhanced_time = static_cast<double>(std::accumulate(
            enhanced_times.begin(), enhanced_times.end(), 0))
                        / static_cast<double>(num_iterations);
    std::cout << "Ignore " << calculated_results[0] << std::endl;
    return ComparisonBenchmarkResult(full_time, reduced_time, enhanced_time);
}

void PrintSurfaceDeboorResult(ComparisonBenchmarkResult& result)
{
    std::cout << "Full : " << result.FirstAlg() << std::endl;
    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
    std::cout << "Enhanced : " << result.ThirdAlg() << std::endl;
    std::cout << "Difference F/E: " << result.Ratio() << std::endl;
}

void PrintCurveDeboorResult(ComparisonBenchmarkResult& result)
{
    std::cout << "Full : " << result.FirstAlg() << std::endl;
    std::cout << "Reduced : " << result.SecondAlg() << std::endl;
    std::cout << "Difference F/R: " << result.Ratio() << std::endl;
}

int main()
{
    bool optimized_tridiagonals = true;
    while (true)
    {
        //std::cout << clock();
        // Console clear ...
        // ... for Windows,
        system("cls");
        // ... for Linux/Unix.
        //system("clear");
        std::cout << "1: Instructions benchmark." << std::endl;
        std::cout << "2: Spline curve benchmark." << std::endl;
        std::cout << "3: Spline surface benchmark." << std::endl;
        std::cout << "4: Spline surface benchmark (in parallel)." << std::endl;
        //std::cout << "5: Compare Csaba T. vs. Vilo K. LU decomposition." <<
        //std::endl;
        //std::cout << "B: Disable/enable optimized LU decomposition in benchmarks."
        //	<< std::endl;
        //std::cout << "B: Enable/disable Vilo Kacala LU decomposition in benchmarks."
        //<< std::endl;
        std::cout << "Q: End program" << std::endl;
        char input;
        std::cin >> input;
        std::cin.get();
        std::cout << std::endl << "---------------" << std::endl;
        unsigned int num_iterations;
        unsigned int num_knots;
        ComparisonBenchmarkResult result(1, 1);
        switch (input)
        {
            case '1':
                std::cout << "Instructions benchmark" << std::endl <<
                          std::endl;
                MulDivBenchmark();
                break;
            case '2':
                std::cout << "Spline curve benchmark" << std::endl << std::endl;
                std::cout << "Enter number of iterations: " << std::endl;
                std::cin >> num_iterations;
                std::cout << "Enter number of knots: " << std::endl;
                std::cin >> num_knots;
                std::cin.get();
                result = CurveBenchmark(num_iterations, num_knots,
                                        optimized_tridiagonals);
                PrintCurveDeboorResult(result);
                break;
            case '3':
                std::cout << "Spline surface benchmark" << std::endl << std::endl;
                std::cout << "Enter number of iterations: " << std::endl;
                std::cin >> num_iterations;
                std::cout << "Enter number of knots: " << std::endl;
                std::cin >> num_knots;
                std::cin.get();

                result = SurfaceBenchmark(num_iterations, num_knots,
                                          false,optimized_tridiagonals);
                PrintSurfaceDeboorResult(result);
                break;
            case '4':
                std::cout << "Parallel spline surface benchmark" << std::endl <<
                          std::endl;
                std::cout << "Enter number of iterations: " << std::endl;
                std::cin >> num_iterations;
                std::cout << "Enter number of knots: " << std::endl;
                std::cin >> num_knots;
                std::cin.get();
                result = SurfaceBenchmark(num_iterations, num_knots, true, optimized_tridiagonals);
                PrintSurfaceDeboorResult(result);
                break;
                /*case '5':
                    LUComparison();
                    break;*/
            case 'q':
            case 'Q':
                return 0;
                /*case 'b':
                case 'B':
                    optimized_tridiagonals = !optimized_tridiagonals;
                    if(optimized_tridiagonals)
                    {
                        std::cout << "Buffering is enabled." << std::endl;
                    }
                    else
                    {
                        std::cout << "Buffering is disabled." << std::endl;
                    }
                    break;*/
        }

        std::cout << "===================" << std::endl;
        /*std::cout << "any key: Restart program." << std::endl;
        std::cout << "Q: End program" << std::endl;*/

        system("pause");
    }
    return 0;
}