#include "stdafx.h"
#include "ComparisonBenchmarkResult.h"


ComparisonBenchmarkResult::ComparisonBenchmarkResult(unsigned long long 
	first_alg_time, unsigned long long second_alg_time)
	:first_alg_(first_alg_time), second_alg_(second_alg_time), ratio_(
		static_cast<double>(first_alg_time) / static_cast<double>(
		 second_alg_time))
{
}

unsigned long long ComparisonBenchmarkResult::FirstAlg() const
{
	return first_alg_;
}

unsigned long long ComparisonBenchmarkResult::SecondAlg() const
{
	return second_alg_;
}

double ComparisonBenchmarkResult::Ratio() const
{
	return ratio_;
}
