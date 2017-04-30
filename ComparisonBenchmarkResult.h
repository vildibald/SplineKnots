#pragma once

class ComparisonBenchmarkResult {
    double first_alg_;
    double second_alg_;
    double ratio_;

public:
    ComparisonBenchmarkResult(double first_alg_time,
                              double second_alg_time);

    unsigned long long FirstAlg() const;

    unsigned long long SecondAlg() const;

    double Ratio() const;
};
