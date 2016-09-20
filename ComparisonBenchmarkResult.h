#pragma once
class ComparisonBenchmarkResult
{
	unsigned long long first_alg_;
	unsigned long long second_alg_;
	double ratio_;

public:
	ComparisonBenchmarkResult(unsigned long long first_alg_time, 
		unsigned long long second_alg_time);
	
	unsigned long long FirstAlg() const;
	
	unsigned long long SecondAlg() const;
	
	double Ratio() const;
};
