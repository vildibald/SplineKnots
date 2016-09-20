#include "stdafx.h"


namespace utils
{
	template <typename Iterator, typename Function>
	void For(Iterator from, Iterator to, Function function, Iterator increment_by, bool in_parallel)
	{
		if (in_parallel)
		{
#pragma omp parallel for
			for (Iterator i = from; i < to; i += increment_by)
			{
				function(i);
			}
		}
		else
		{
			// Loop for sequential computations. 
			// '#pragma omp parallel for if(in_parallel)' in above loop will execute that loop in one thread, 
			// but still with overhead of OpenMP thread creation. 
			for (Iterator i = from; i < to; i += increment_by)
			{
				function(i);
			}
		}
	}

	template <typename T>
	T* InitArray(size_t length, T arrayToInit[], T value)
	{
		for (size_t i = 0; i < length; i++)
		{
			arrayToInit[i] = value;
		}
		return arrayToInit;
	}

	template <typename T>
	void DeleteJaggedArray(T** jaggedArray, size_t rows, size_t columns)
	{
		for (size_t i = 0; i < rows; i++)
		{
			delete[] jaggedArray[i];
			jaggedArray[i] = nullptr;
		}
		delete[] jaggedArray;
		jaggedArray = nullptr;
	}

	template <typename T>
	T** CreateJaggedArray(size_t rows, size_t columns)
	{
		auto res = new T*[rows];
		for (size_t i = 0; i < columns; i++)
		{
			res[i] = new T[columns];
		}

		return res;
	}

	template <typename T>
	double Average(T arr[], size_t arr_size)
	{
		T sum = 0;
		for (size_t i = 0; i < arr_size; i++)
		{
			sum += arr[i];
		}
		return static_cast<double>(sum) / arr_size;
	}
}
