#pragma once
#include <omp.h>
#include <thread>
#include <vector>

namespace utils
{
	extern unsigned int num_threads;

	template <typename Iterator, typename Function>
	void For(Iterator from, Iterator to, Function function,
		Iterator increment_by = 1, bool in_parallel = false);

	template <typename T>
	T* InitArray(size_t length, T arrayToInit[], T value);
	
	template <typename T>
	void DeleteJaggedArray(T** jaggedArray, size_t rows, size_t columns);
	
	template <typename T>
	T** CreateJaggedArray(size_t rows, size_t columns);
	
	std::vector<double> SolveCsabaDeboorTridiagonalSystem(double b,
		double right_side[], unsigned int num_equations,
		double last_main_diagonal_value = DBL_MIN);
	
	void SolveCsabaDeboorTridiagonalSystemBuffered(double main_diagonal_value,
		double right_side[], unsigned int num_equations, double L[],
		double U[], double Y[], double D[],
		double last_main_diagonal_value = DBL_MIN);
	
	void SolveTridiagonalSystem(double lower_diagonal[],
		double main_diagonal[], double upper_diagonal[], double right_side[],
		size_t num_equations);
	
	void SolveTridiagonalSystemBuffered(double lower_diagonal[],
		double main_diagonal[],
		double upper_diagonal[], double right_side[], size_t num_equations,
		double buffer[]);
	
	void SolveDeboorTridiagonalSystem(double lower_diagonal_value,
		double main_diagonal_value, double upper_diagonal_value,
		double right_side[], size_t num_equations,
		double last_main_diagonal_value = DBL_MIN);
	
	void SolveDeboorTridiagonalSystemBuffered(double lower_diagonal_value,
		double main_diagonal_value, double upper_diagonal_value,
		double right_side[], size_t num_equations, double buffer[],
		double last_main_diagonal_value = DBL_MIN);
	
	void SolveDeboorTridiagonalSystem(double main_diagonal_value,
		double right_side[], size_t num_equations,
		double last_main_diagonal_value = DBL_MIN);
	
	void SolveDeboorTridiagonalSystemBuffered(double main_diagonal_value,
		double right_side[], size_t num_equations, double buffer[], double
		last_main_diagonal_value = DBL_MIN);
	
	template <typename T>
	double Average(T arr[], size_t arr_size);
};

#include "utils_template.cpp"