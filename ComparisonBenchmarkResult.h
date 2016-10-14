#pragma once

class ComparisonBenchmarkResult {
    unsigned long long first_alg_;
    unsigned long long second_alg_;
    unsigned long long third_alg_;
    double ratio_;

public:
    ComparisonBenchmarkResult(unsigned long long first_alg_time,
                              unsigned long long second_alg_time,
                              unsigned long long third_alg_time = 0);

    unsigned long long FirstAlg() const;

    unsigned long long SecondAlg() const;

    unsigned long long ThirdAlg() const;

    double Ratio() const;
};
