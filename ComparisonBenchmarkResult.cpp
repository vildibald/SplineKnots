#include "stdafx.h"
#include "ComparisonBenchmarkResult.h"


ComparisonBenchmarkResult::ComparisonBenchmarkResult(unsigned long long
                                                     firstAlgTime,
                                                     unsigned long long secondAlgTime,
                                                     unsigned long long
                                                     thirdAlgTime)
        : firstAlg(firstAlgTime), secondAlg(secondAlgTime), thirdAlg(
        thirdAlgTime) {
    auto divisor = thirdAlgTime == -1 ? secondAlgTime : thirdAlgTime;
    ratio = static_cast<double>(firstAlgTime) / static_cast<double>(
            divisor);
}

unsigned long long ComparisonBenchmarkResult::FirstAlg() const {
    return firstAlg;
}

unsigned long long ComparisonBenchmarkResult::SecondAlg() const {
    return secondAlg;
}

unsigned long long ComparisonBenchmarkResult::ThirdAlg() const {
    return thirdAlg;
}

double ComparisonBenchmarkResult::Ratio() const {
    return ratio;
}
