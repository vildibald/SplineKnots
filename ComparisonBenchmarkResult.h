#pragma once

class ComparisonBenchmarkResult {
    unsigned long long firstAlg;
    unsigned long long secondAlg;
    unsigned long long thirdAlg;
    double ratio;

public:
    ComparisonBenchmarkResult(unsigned long long firstAlgTime,
                              unsigned long long secondAlgTime,
                              unsigned long long thirdAlgTime = 0);

    unsigned long long FirstAlg() const;

    unsigned long long SecondAlg() const;

    unsigned long long ThirdAlg() const;

    double Ratio() const;
};
