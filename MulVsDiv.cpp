#include "stdafx.h"
#include "MulVsDiv.h"
#include <ostream>
#include <iostream>
#include <locale>
#include <vector>
#include <numeric>
#include <windows.h>
#include <random>
#include "StopWatch.h"

void 
MulVsDiv::ResetArrays(const int length, double* a, double* b, double& ignoreit)
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	ignoreit += a[((int)dist(mt) % (int)(length))] + b[((int)dist(mt) % (int)(length))];
	for(size_t i = 0; i < length; i++)
	{
		a[i] = 6 * (dist(mt) / (to)) + 2;
		b[i] = 2 * (dist(mt) / (to)) + 1;
	}
	ignoreit += a[(int)dist(mt) % (int)(length)] + b[(int)dist(mt) % (int)(length)];
	ignoreit = 1 / ignoreit;
}

void
MulVsDiv::ResetMatrix(const int rows, const int columns, double** matrix, 
	double& ignoreit)
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	ignoreit += matrix[(int)(dist(mt)) % (int)(rows)][(int)dist(mt) % (int)(columns)];

	for(size_t i = 0; i < rows; i++)
	{
		for(size_t j = 0; j < columns; j++)
		{
			matrix[i][j] = 2 * ((int)(dist(mt)) / (to)) + 1;
		}
	}
	ignoreit += matrix[(int)(dist(mt)) % (int)(rows)][(int)dist(mt) % (int)(columns)];
	ignoreit = 1 / ignoreit;
}

void
MulVsDiv::ResetList(const int length, std::forward_list<double>& list, double& ignoreit)
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	list.clear();
	for(size_t i = 0; i < length; i++)
	{
		list.push_front(6 * ((int)(dist(mt)) / (to)) + 2);
	}
	ignoreit += list.front();
	ignoreit = 1 / ignoreit;
}

void
MulVsDiv::Loop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Simple loop:\n---" << std::endl;
	double a[length], b[length], c[length];
	//double d[length], e[length], f[length], g[length];
	//double d1[length], e1[length], f1[length], g1[length],x[length];
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);
	/*ResetArrays(length, d, e, ignoreit);
	ResetArrays(length, f, g, ignoreit);
	ResetArrays(length, g, x, ignoreit);*/
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will 
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and 
		// nonvectorized comparison i specifically disabled vectorization in 
		// this method
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] + b[i];// +d[i] + e[i] +f[i] + g[i] + d1[i] + e1[i] + f1[i] + g1[i] + x[i];
		}
	}
	sw.Stop();
	auto add_time = sw.EllapsedTime();
	std::cout << "Addition: " << add_time << std::endl;
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, b, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			a[i] = c[i] * b[i];
		}
	}
	sw.Stop();
	auto mul_time = sw.EllapsedTime();
	std::cout << "Multiplication: " << mul_time << std::endl;
	ignoreit -= a[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, b, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			b[i] = a[i] / c[i];
		}
	}
	sw.Stop();
	auto div_time = sw.EllapsedTime();
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];
	std::cout << "Division: " << div_time << std::endl;
	ignoreit -= a[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, b, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = 1 / a[i];
		}
	}
	sw.Stop();
	auto recdiv_time = sw.EllapsedTime();
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];
	std::cout << "Reciprocal division: " << recdiv_time << std::endl;
	/*std::cout << "Addition faster than multiplication: " << static_cast<double>
		(mul_time) / static_cast<double>(add_time) << std::endl;
	std::cout << "Multiplication faster than division: " << static_cast<double>
		(div_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Multiplication faster than reciprocal division: " << static_cast<double>
		(recdiv_time) / static_cast<double>(mul_time) << std::endl;*/
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::ArrayAndNumberLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Simple array and number loop:\n---" << std::endl;
	double a[length], b, c[length];
	auto ignoreit = 0.0;
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	ResetArrays(length, a, c, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] + b;
		}
	}
	sw.Stop();
	auto add_time = sw.EllapsedTime();
	std::cout << "Addition: " << add_time << std::endl;
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, c, ignoreit);
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			a[i] = c[i] * b;
		}
	}
	sw.Stop();
	auto mul_time = sw.EllapsedTime();
	std::cout << "Multiplication: " << mul_time << std::endl;
	ignoreit -= a[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, c, ignoreit);
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] / b;
		}
	}
	sw.Stop();
	auto div_time = sw.EllapsedTime();
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b;
	std::cout << "Division: " << div_time << std::endl;
	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::LoopVectorized()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 256;
	const int loops = 1e6;
	std::cout << "Vectorized loop:\n---" << std::endl;
	double a[length], b[length], c[length];
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);
	int l = 0;
	auto start = clock();
add:
	for(int i = 0; i < length; i++)
	{
		b[i] = a[i] + c[i];
	}
	// MSVC doesn't vectorize nested loops (message 1300 - too little
	// computation to vectorize) as mentioned on line 23 in function 'Loop'.
	// However if nested loop is replaced with this nasty workaround SIMD
	// vectorization will happen. Same condition apply for mul/div/rcp loops
	while(l < loops)
	{
		++l;
		goto add;
	}
	auto add_time = clock() - start;
	std::cout << "Addition: " << add_time << std::endl;
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	ResetArrays(length, a, b, ignoreit);

	l = 0;
	start = clock();
mul:
	for(int i = 0; i < length; ++i)
	{
		c[i] = a[i] * b[i];
	}
	while(l < loops)
	{
		++l;
		goto mul;
	}
	auto mul_time = clock() - start;
	std::cout << "Multiplication: " << mul_time << std::endl;
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	ResetArrays(length, a, b, ignoreit);

	l = 0;
	start = clock();
div:
	for(int i = 0; i < length; ++i)
	{
		a[i] = c[i] / b[i];
	}
	while(l < loops)
	{
		++l;
		goto div;
	}
	auto div_time = clock() - start;
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	std::cout << "Division: " << div_time << std::endl;
	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::DynamicArrayLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Dynamic array loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length), cv(length);
	//std::vector<double> dv(length), ev(length), fv(length), gv(length);
	//std::vector<double> d1v(length), e1v(length), f1v(length), g1v(length), xv(length);
	double *a = &av.front(), *b = &bv.front(), *c = &cv.front();
	//double *d = &dv.front(), *e = &ev.front(), *f = &fv.front(), *g = &gv.front();
	//double *d1 = &d1v.front(), *e1 = &e1v.front(), *f1 = &f1v.front(), *g1 = &g1v.front(),*x=&xv.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);
	//ResetArrays(length, d, e, ignoreit);
	//ResetArrays(length, f, g, ignoreit);
	//ResetArrays(length, g, x, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC doesn't vectorize nested loops (message 1300 - too little
		// computation to vectorize) as mentioned on line 23 in function 'Loop'.
		// However if nested loop is replaced with this nasty workaround SIMD
		// vectorization will happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] + b[i];// *d[i] + (b[i] + d[i])*a[i];// +d[i] + e[i] + f[i] + g[i] + d1[i] + e1[i] + f1[i] + g1[i] + x[i];
		}
	}
	sw.Stop();
	auto add_time = sw.EllapsedTime();
	std::cout << "Addition: " << add_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] * b[i];// *d[i] * e[i] * f[i] * g[i] * d1[i] * e1[i] * f1[i] * g1[i] * x[i];
		}
	}
	sw.Stop();
	auto mul_time = sw.EllapsedTime();
	std::cout << "Multiplication: " << mul_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = b[i] / a[i];// / d[i] / e[i] / f[i] / g[i] / d1[i] / e1[i] / f1[i] / g1[i] / x[i];
		}
	}
	sw.Stop();
	auto div_time = sw.EllapsedTime();
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	std::cout << "Division: " << div_time << std::endl;
	ignoreit -= a[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a, b, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = 1 / a[i];
		}
	}
	sw.Stop();
	auto recdiv_time = sw.EllapsedTime();
	ignoreit -= c[((int)(dist(mt)) % static_cast<int>(length))];
	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];
	std::cout << "Reciprocal division: " << recdiv_time << std::endl;
	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	//std::cout << "Multiplication faster than reciprocal division: " << static_cast<double>
	//	(recdiv_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::DynamicListLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Dynamic list loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length);
	double *a = &av.front(), *b = &bv.front();
	std::forward_list<double> cl;
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);
	ResetList(length, cl, ignoreit);


	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC doesn't vectorize nested loops (message 1300 - too little
		// computation to vectorize) as mentioned on line 23 in function 'Loop'.
		// However if nested loop is replaced with this nasty workaround SIMD
		// vectorization will happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		int i = 0;
		for(auto& c : cl)
		{
			c = a[i] + b[i];
			++i;
		}
	}
	sw.Stop();
	auto add_time = sw.EllapsedTime();
	std::cout << "Addition: " << add_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	auto it = cl.begin();
	std::advance(it, ((int)(dist(mt)) % (int)(length)));
	ignoreit /= *it;

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		int i = 0;
		for(auto& c : cl)
		{
			c = a[i] * b[i];
			++i;
		}
	}
	sw.Stop();
	auto mul_time = sw.EllapsedTime();
	std::cout << "Multiplication: " << mul_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	it = cl.begin();
	std::advance(it, ((int)(dist(mt)) % (int)(length)));
	ignoreit /= *it;
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		int i = 0;
		for(auto& c : cl)
		{
			c = b[i] / a[i];
			++i;
		}
	}
	sw.Stop();
	auto div_time = sw.EllapsedTime();
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	it = cl.begin();
	std::advance(it, ((int)(dist(mt)) % (int)(length)));
	ignoreit /= *it;
	std::cout << "Division: " << div_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	it = cl.begin();
	std::advance(it, ((int)(dist(mt)) % (int)(length)));
	ignoreit /= *it;
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		int i = 0;
		for(auto& c : cl)
		{
			c = 1.0 / a[i];
			++i;
		}
	}
	sw.Stop();
	auto recdiv_time = sw.EllapsedTime();
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	it = cl.begin();
	std::advance(it, ((int)(dist(mt)) % (int)(length)));
	ignoreit /= *it;
	std::cout << "Reciprocal division: " << recdiv_time << std::endl;
	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	//std::cout << "Multiplication faster than reciprocal division: " << static_cast<double>
	//	(recdiv_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::DynamicArrayAndNumberLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Dynamic array with number loop:\n---" << std::endl;
	std::vector<double> av(length), cv(length);
	double *a = &av.front(), *c = &cv.front();
	double b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;

	auto ignoreit = 0.0;

	ResetArrays(length, a, c, ignoreit);

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC doesn't vectorize nested loops (message 1300 - too little
		// computation to vectorize) as mentioned on line 23 in function 'Loop'.
		// However if nested loop is replaced with this nasty workaround SIMD
		// vectorization will happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] + b;// +d[i] + e[i] + f[i] + g[i] + d1[i] + e1[i] + f1[i] + g1[i] + x[i];
		}
	}
	sw.Stop();
	auto add_time = sw.EllapsedTime();
	std::cout << "Addition: " << add_time << std::endl;

	ResetArrays(length, a, c, ignoreit);
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] * b;// *d[i] * e[i] * f[i] * g[i] * d1[i] * e1[i] * f1[i] * g1[i] * x[i];
		}
	}
	sw.Stop();
	auto mul_time = sw.EllapsedTime();
	std::cout << "Multiplication: " << mul_time << std::endl;
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	ResetArrays(length, a, c, ignoreit);
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = b / a[i];// / d[i] / e[i] / f[i] / g[i] / d1[i] / e1[i] / f1[i] / g1[i] / x[i];
		}
	}
	sw.Stop();
	auto div_time = sw.EllapsedTime();
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b;
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	std::cout << "Division number/vector: " << div_time << std::endl;
	b = 6 * ((double)(int)(dist(mt)) / (to)) + 2;
	ResetArrays(length, a, c, ignoreit);
	ignoreit -= c[((int)(dist(mt)) % (int)(length))];
	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			c[i] = a[i] / b;// / d[i] / e[i] / f[i] / g[i] / d1[i] / e1[i] / f1[i] / g1[i] / x[i];
		}
	}
	sw.Stop();
	auto div_time1 = sw.EllapsedTime();
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b;
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	std::cout << "Division: " << div_time << std::endl;
	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than num/vector division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	//std::cout << "Multiplication faster than vector/num division: " << static_cast<double>
	//	(div_time1) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::DynamicArrayLoopVectorized()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 256;
	const int loops = 1e6;
	std::cout << "Vectorized dynamic array loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length), cv(length);
	double *a = &av.front(), *b = &bv.front(), *c = &bv.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);

	int l = 0;
	auto start = clock();
add:
	for(int i = 0; i < length; i++)
	{
		c[i] = a[i] + b[i];
	}
	// MSVC doesn't vectorize nested loops (message 1300 - too little
	// computation to vectorize) in function 'DynamicArrayLoop'. However if
	// nested loop is replaced with this nasty workaround SIMD vectorization
	// will happen. Same condition apply for mul/div/rcp loops
	while(l < loops)
	{
		++l;
		goto add;
	}
	auto add_time = clock() - start;
	std::cout << "Addition: " << add_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	l = 0;
	start = clock();
mul:
	for(int i = 0; i < length; ++i)
	{
		c[i] = a[i] * b[i];
	}
	while(l < loops)
	{
		++l;
		goto mul;
	}
	auto mul_time = clock() - start;
	std::cout << "Multiplication: " << mul_time << std::endl;

	ResetArrays(length, a, b, ignoreit);
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	l = 0;
	start = clock();
div:
	for(int i = 0; i < length; ++i)
	{
		c[i] = a[i] / b[i];
	}
	while(l < loops)
	{
		++l;
		goto div;
	}
	auto div_time = clock() - start;
	ignoreit /= c[((int)(dist(mt)) % (int)(length))];
	ignoreit += a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	std::cout << "Division: " << div_time << std::endl;

	//std::cout << "Addition faster than multiplication: " << static_cast<double>
	//	(mul_time) / static_cast<double>(add_time) << std::endl;
	//std::cout << "Multiplication faster than division: " << static_cast<double>
	//	(div_time) / static_cast<double>(mul_time) << std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::CsabaDynamicArrayLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 256;//1024 *1024* 8;
	const int loops = 1e6;//1e3 / 2;
	std::cout << "Csaba loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length);
	double *a = &av.front(), *b = &bv.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);
	auto idx1 = ((int)(dist(mt)) % static_cast<int>(length)), idx2 = ((int)(dist(mt)) %
		static_cast<int>(length - 1));
	auto start = clock();
	for(size_t l = 0; l < loops; l++)
	{
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		// happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparisoni specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for(int i = 1; i < length - 1; ++i)
		{
			a[i] = b[idx1] + b[idx2];
		}
	}
	auto add_time = clock() - start;
	std::cout << "Addition: " << add_time << std::endl;

	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];

	idx1 = ((int)(dist(mt)) % static_cast<int>(length)) , idx2 = ((int)(dist(mt)) %
		static_cast<int>(length));

	start = clock();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length - 1; ++i)
		{
			a[i] = b[idx1] * b[idx2];
		}
	}
	auto mul_time = clock() - start;
	std::cout << "Multiplication: " << mul_time << std::endl;

	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];
	idx1 = ((int)(dist(mt)) % static_cast<int>(length)) , idx2 = ((int)(dist(mt)) %
		static_cast<int>(length));
	start = clock();
	for(size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for(int i = 0; i < length; ++i)
		{
			a[i] = b[idx1] / b[idx2];
		}
	}
	auto div_time = clock() - start;
	ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] + b[((int)(dist(mt)) %
		static_cast<int>(length))];
	std::cout << "Division: " << div_time << std::endl;
	//std::cout << "Addition faster than multiplication: " <<
	//	static_cast<double>(mul_time) / static_cast<double>(add_time) <<
	//	std::endl;
	//std::cout << "Multiplication faster than division: " <<
	//	static_cast<double>(div_time) / static_cast<double>(mul_time) <<
	//	std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::BackwardDependendDynamicArrayLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Backward dependend dynamic loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length);
	double *a = &av.front(), *b = &bv.front();

	std::vector<long> add_times;
	std::vector<long> mul_times;
	std::vector<long> div_times;
	add_times.reserve(loops);
	mul_times.reserve(loops);
	div_times.reserve(loops);
	auto ignoreit = 0.0;

	for(size_t i = 0; i < length; i++)
	{
		av[i] = 8 * (static_cast<double>((int)(dist(mt))) / (to)) + 1;
		bv[i] = DBL_MIN ;
	}

	//Division

	for(size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for(size_t j = 1; j < length; j++)
		{
			a[j] = a[j] / a[j - 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
		//ignoreit[0] += std::accumulate(bv.begin(), bv.end(), 0);
	}

	auto div_time = static_cast<double>(std::accumulate(div_times.begin(),
		div_times.end(), 0));
	std::cout << "Division:\t\t" << div_time << std::endl;
	//Multiplication

	for(size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for(size_t j = 1; j < length; j++)
		{
			a[j] = a[j] * a[j - 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto mul_time = static_cast<double>(std::accumulate(mul_times.begin(),
		mul_times.end(), 0));
	std::cout << "Multiplication:\t\t" << mul_time << std::endl;
	//Addition

	for(size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for(size_t j = 1; j < length; j++)
		{
			a[j] = a[j] + a[j - 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto add_time = static_cast<double>(std::accumulate(add_times.begin(),
		add_times.end(), 0));
	std::cout << "Addition:\t\t" << add_time << std::endl;

	//std::cout << "Addition faster than multiplication:\t" <<
	//	static_cast<double>(mul_time) / static_cast<double>(add_time) <<
	//	std::endl;
	//std::cout << "Multiplication faster than division:\t" <<
	//	static_cast<double>(div_time) / static_cast<double>(mul_time) <<
	//	std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::ForwardDependendDynamicArrayLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Forward dependend dynamic loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length);
	double *a = &av.front(), *b = &bv.front();

	std::vector<long> add_times;
	std::vector<long> mul_times;
	std::vector<long> div_times;
	add_times.reserve(loops);
	mul_times.reserve(loops);
	div_times.reserve(loops);
	auto ignoreit = 0.0;

	for(size_t i = 0; i < length; i++)
	{
		av[i] = 8 * (static_cast<double>((int)(dist(mt))) / (to)) + 1;
		bv[i] = DBL_MIN ;
	}

	//Division

	for(size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for(size_t j = 0; j < length-1; j++)
		{
			a[j] = a[j] / a[j + 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
		//ignoreit[0] += std::accumulate(bv.begin(), bv.end(), 0);
	}

	auto div_time = static_cast<double>(std::accumulate(div_times.begin(),
		div_times.end(), 0));
	std::cout << "Division:\t\t" << div_time << std::endl;
	//Multiplication

	for(size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++)
		{
			a[j] = a[j] * a[j + 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto mul_time = static_cast<double>(std::accumulate(mul_times.begin(),
		mul_times.end(), 0));
	std::cout << "Multiplication:\t\t" << mul_time << std::endl;
	//Addition

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++)
		{
			a[j] = a[j]+ a[j + 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto add_time = static_cast<double>(std::accumulate(add_times.begin(),
		add_times.end(), 0));
	std::cout << "Addition:\t\t" << add_time << std::endl;

	std::cout << "Addition faster than multiplication:\t" <<
		static_cast<double>(mul_time) / static_cast<double>(add_time) <<
		std::endl;
	std::cout << "Multiplication faster than division:\t" <<
		static_cast<double>(div_time) / static_cast<double>(mul_time) <<
		std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void
MulVsDiv::DependendDynamicArrayLoop()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 512;
	const int loops = 1e6 / 4;
	std::cout << "Dependend dynamic loop:\n---" << std::endl;
	std::vector<double> av(length);
	double *a = &av.front();

	std::vector<long> add_times;
	std::vector<long> mul_times;
	std::vector<long> div_times;
	add_times.reserve(2*loops);
	mul_times.reserve(2*loops);
	div_times.reserve(2*loops);

	auto ignoreit = 0.0;

	for (size_t i = 0; i < length; i++)
	{
		av[i] = 8 * (static_cast<double>((int)(dist(mt))) / (to)) + 1;
	}
	//  Backward  ///////////////////////////////////////////////////////////////////////

	//Division

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 1; j < length; j++)
		{
			a[j] = a[j] / a[j - 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 16;
		//ignoreit[0] += std::accumulate(bv.begin(), bv.end(), 0);
	}

	//Multiplication

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 1; j < length; j++)
		{
			a[j] = a[j] * a[j - 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 16;
	}

	//Addition

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 1; j < length; j++)
		{
			a[j] = a[j] + a[j - 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 16;
	}

	//  Forward  ///////////////////////////////////////////////////////////////////////

	for (size_t i = 0; i < length; i++)
	{
		av[i] = 8 * (static_cast<double>((int)(dist(mt))) / (to)) + 1;

	}

	//Division

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++)
		{
			a[j] = a[j] / a[j + 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto div_time = static_cast<double>(std::accumulate(div_times.begin(),
		div_times.end(), 0));
	std::cout << "Division:\t\t" << div_time << std::endl;
	//Multiplication

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++)
		{
			a[j] = a[j] * a[j + 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto mul_time = static_cast<double>(std::accumulate(mul_times.begin(),
		mul_times.end(), 0));
	std::cout << "Multiplication:\t\t" << mul_time << std::endl;
	//Addition

	for (size_t i = 0; i < loops; i++)
	{
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++)
		{
			a[j] = a[j] + a[j + 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int)(dist(mt)) % static_cast<int>(length))] / 8;
	}

	auto add_time = static_cast<double>(std::accumulate(add_times.begin(),
		add_times.end(), 0));
	std::cout << "Addition:\t\t" << add_time << std::endl;

	std::cout << "Addition faster than multiplication:\t" <<
		static_cast<double>(mul_time) / static_cast<double>(add_time) <<
		std::endl;
	std::cout << "Multiplication faster than division:\t" <<
		static_cast<double>(div_time) / static_cast<double>(mul_time) <<
		std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}


void
MulVsDiv::MemCpy()
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
#pragma optimize("", off)
	StopWatch sw;
	const int length = 512;
	const int loops = 1e6 / 2;
	std::cout << "Memory copy loop:\n---" << std::endl;
	std::vector<double> av(length), bv(length), cv(length * 2);
	double *a = &av.front(), *b = &bv.front(), *c = &cv.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a, b, ignoreit);

	ignoreit /= a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	sw.Start();
	//#pragma optimize("", off)
	for(size_t l = 0; l < loops; l++)
	{
		memcpy(a, b, length);
	}
	//#pragma optimize("", on)
	sw.Stop();
	auto con_time = sw.EllapsedTime();
	std::cout << "Continuous:\t\t" << con_time << std::endl;
	ignoreit /= a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	ResetArrays(length, a, b, ignoreit);
	ignoreit /= a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];

	sw.Start();
	for(size_t l = 0; l < loops; l++)
	{
		//#pragma optimize("", off)
		for(int i = 0, k = 0; i < length; ++i , k += 2)
		{
			c[k] = a[i];
		}
		//#pragma optimize("", on)
	}
	sw.Stop();
	ignoreit /= a[((int)(dist(mt)) % (int)(length))] + b[((int)(dist(mt)) % (int)(length))];
	auto uncon_time = sw.EllapsedTime();
	std::cout << "Uncontinuous:\t\t" << uncon_time << std::endl;


	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
#pragma optimize("", on)
}

void MulVsDiv::ArithmeticInstructionParallelism()
{
	std::cout << "Instruction parallelism test:\n---" << std::endl;
	double add1_time, mul1_time, div1_time;
	ArithmeticInstructionParallelismSingleOperand(&add1_time, &mul1_time,
		&div1_time);
	double add10_time, mul10_time, div10_time;
	ArithmeticInstructionParallelismTenOperands(&add10_time, &mul10_time,
		&div10_time);
	std::cout << "Times with 1 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add1_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul1_time << std::endl;
	std::cout << "Division:\t\t" << div1_time <<"\n"<< std::endl;

	std::cout << "Times with 10 operations per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add10_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul10_time << std::endl;
	std::cout << "Division:\t\t" << div10_time << std::endl;

}

void MulVsDiv::ArithmeticInstructionParallelismSingleOperand(double* add_time,
	double* mul_time, double* div_time)
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length);
	double *a0 = &a0v.front(), *a1 = &a1v.front(), *a2 = &a2v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a1, a2, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] + a2[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a1, a2, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] * a2[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a1, a2, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] / a2[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
}


void MulVsDiv::ArithmeticInstructionParallelismTenOperands(double* add_time,
	double* mul_time, double* div_time)
{
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
		a5v(length), a6v(length), a7v(length), a8v(length), a9v(length),
		a10v(length), a11v(length);
	double *a0 = &a0v.front(), *a1 = &a1v.front(), *a2 = &a2v.front(), *a3 = &a3v.front(), *a4 = &a4v.front(),
		*a5 = &a5v.front(), *a6 = &a6v.front(), *a7 = &a7v.front(), *a8 = &a8v.front(), *a9 = &a9v.front(),
		*a10 = &a10v.front(), *a11 = &a11v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
				a8[i] + a9[i] + a10[i] + a11[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] * a2[i] * a3[i] * a4[i] * a5[i] * a6[i] * a7[i] *
				a8[i] * a9[i] * a10[i] * a11[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++)
	{
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i)
		{
			a0[i] = a1[i] / a2[i] / a3[i] / a4[i] / a5[i] / a6[i] / a7[i]
				/ a8[i] / a9[i] / a10[i] / a11[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int)(dist(mt)) % static_cast<int>(length))];
}

void
MulVsDiv::BenchAll()
{
	//Loop();
	//ArrayAndNumberLoop();
	//LoopVectorized();
	//DynamicArrayLoop();	
	//DynamicArrayAndNumberLoop();
	//DynamicListLoop();
	MemCpy();
	//DynamicArrayLoopVectorized();
	//ForwardDependendDynamicArrayLoop();
	//BackwardDependendDynamicArrayLoop();
	DependendDynamicArrayLoop();
	ArithmeticInstructionParallelism();
	//CsabaDynamicArrayLoop();
}

MulVsDiv::MulVsDiv()
{
}


MulVsDiv::~MulVsDiv()
{
}
