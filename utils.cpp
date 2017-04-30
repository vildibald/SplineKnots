#include "stdafx.h"
#include "utils.h"
#include <omp.h>
#include <vector>

namespace utils
{
	unsigned int num_threads = std::thread::hardware_concurrency();

	void SolveDeboorTridiagonalSystemBuffered(double main_diagonal_value,
		double right_side[], size_t num_equations, double buffer[], double
		last_main_diagonal_value)
	{
		if (last_main_diagonal_value == DBL_MIN)
			last_main_diagonal_value = main_diagonal_value;
		auto m0 = 1 / main_diagonal_value;
		buffer[0] = m0;
		right_side[0] *= m0;
		auto lastindex = num_equations - 1;
		for (size_t i = 1; i < lastindex; i++)
		{
			auto m = 1 / (main_diagonal_value - buffer[i - 1]);
			buffer[i] = m;
			right_side[i] = (right_side[i] - right_side[i - 1]) * m;
		}

		m0 = 1 / (last_main_diagonal_value - buffer[lastindex - 1]);
		buffer[lastindex] = m0;
		right_side[lastindex] = (right_side[lastindex]
			- right_side[lastindex - 1]) * m0;

		for (size_t i = num_equations - 1; i-- > 0;)
		{
			right_side[i] -= buffer[i] * right_side[i + 1];
		}
	}

	void SolveDeboorTridiagonalSystem(double main_diagonal_value,
		double right_side[], size_t num_equations,
		double last_main_diagonal_value)
	{
		auto buffer = new double[num_equations];
		SolveDeboorTridiagonalSystemBuffered(main_diagonal_value, right_side,
			num_equations, buffer, last_main_diagonal_value);
		delete[] buffer;
	}


}
