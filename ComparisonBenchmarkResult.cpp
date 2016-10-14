#include "stdafx.h"
#include "ComparisonBenchmarkResult.h"


ComparisonBenchmarkResult::ComparisonBenchmarkResult(unsigned long long
                                                     first_alg_time,
                                                     unsigned long long second_alg_time,
                                                     unsigned long long
                                                     third_alg_time)
        : first_alg_(first_alg_time), second_alg_(second_alg_time), third_alg_(
        third_alg_time) {
    auto divisor = third_alg_time == -1 ? second_alg_time : third_alg_time;
    ratio_ = static_cast<double>(first_alg_time) / static_cast<double>(
            divisor);
}

unsigned long long ComparisonBenchmarkResult::FirstAlg() const {
    return first_alg_;
}

unsigned long long ComparisonBenchmarkResult::SecondAlg() const {
    return second_alg_;
}

unsigned long long ComparisonBenchmarkResult::ThirdAlg() const {
    return third_alg_;
}

double ComparisonBenchmarkResult::Ratio() const {
    return ratio_;
}
