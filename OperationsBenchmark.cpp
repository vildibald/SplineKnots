#include "stdafx.h"
#include "OperationsBenchmark.h"
#include <ostream>
#include <iostream>
#include <locale>
#include <vector>
#include <numeric>
#include <windows.h>
#include <random>
#include "StopWatch.h"

void
OperationsBenchmark::ResetArray(const int length, double* a, double& ignoreit) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	ignoreit += a[((int) dist(mt)%(int) (length))];
	for (size_t i = 0; i < length; i++) {
		a[i] = 6*(dist(mt)/(to)) + 2;
	}
	ignoreit += a[(int) dist(mt)%(int) (length)];
	ignoreit = 1/ignoreit;
}

void
OperationsBenchmark::ResetArrays(const int length, double* a, double* b, double& ignoreit) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	ignoreit += a[((int) dist(mt)%(int) (length))] + b[((int) dist(mt)%(int) (length))];
	for (size_t i = 0; i < length; i++) {
		a[i] = 6*(dist(mt)/(to)) + 2;
		b[i] = 2*(dist(mt)/(to)) + 1;
	}
	ignoreit += a[(int) dist(mt)%(int) (length)] + b[(int) dist(mt)%(int) (length)];
	ignoreit = 1/ignoreit;
}

void
OperationsBenchmark::ArithemticOperations() {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	const int length = 512;
	const int loops = 1e6/4;
	std::cout << "Dependend dynamic loop:\n---" << std::endl;
	std::vector<double> av(length);
	double* a = &av.front();

	std::vector<long> add_times;
	std::vector<long> mul_times;
	std::vector<long> div_times;
	add_times.reserve(2*loops);
	mul_times.reserve(2*loops);
	div_times.reserve(2*loops);

	auto ignoreit = 0.0;

	for (size_t i = 0; i < length; i++) {
		av[i] = 8*(static_cast<double>((int) (dist(mt)))/(to)) + 1;
	}
	//  Backward  ///////////////////////////////////////////////////////////////////////

	//Division

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 1; j < length; j++) {
			a[j] = a[j]/a[j - 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/16;
		//ignoreit[0] += std::accumulate(bv.begin(), bv.end(), 0);
	}

	//Multiplication

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 1; j < length; j++) {
			a[j] = a[j]*a[j - 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/16;
	}

	//Addition

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 1; j < length; j++) {
			a[j] = a[j] + a[j - 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/16;
	}

	//  Forward  ///////////////////////////////////////////////////////////////////////

	for (size_t i = 0; i < length; i++) {
		av[i] = 8*(static_cast<double>((int) (dist(mt)))/(to)) + 1;

	}

	//Division

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++) {
			a[j] = a[j]/a[j + 1];
		}
		auto finish = clock();
		div_times.push_back(finish - start);

		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/8;
	}

	auto div_time = static_cast<double>(std::accumulate(div_times.begin(),
														div_times.end(), 0));
	std::cout << "Division:\t\t" << div_time << std::endl;
	//Multiplication

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++) {
			a[j] = a[j]*a[j + 1];
		}
		auto finish = clock();
		mul_times.push_back(finish - start);
		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/8;
	}

	auto mul_time = static_cast<double>(std::accumulate(mul_times.begin(),
														mul_times.end(), 0));
	std::cout << "Multiplication:\t\t" << mul_time << std::endl;
	//Addition

	for (size_t i = 0; i < loops; i++) {
		auto start = clock();
		for (size_t j = 0; j < length - 1; j++) {
			a[j] = a[j] + a[j + 1];
		}
		auto finish = clock();
		add_times.push_back(finish - start);
		ignoreit += a[((int) (dist(mt))%static_cast<int>(length))]/8;
	}

	auto add_time = static_cast<double>(std::accumulate(add_times.begin(),
														add_times.end(), 0));
	std::cout << "Addition:\t\t" << add_time << std::endl;

	std::cout << "Addition faster than multiplication:\t" <<
			  static_cast<double>(mul_time)/static_cast<double>(add_time) <<
			  std::endl;
	std::cout << "Multiplication faster than division:\t" <<
			  static_cast<double>(div_time)/static_cast<double>(mul_time) <<
			  std::endl;
	std::cout << "Just ignore it: " << ignoreit << std::endl << std::endl;
}

void OperationsBenchmark::ArithmeticInstructionParallelism() {
	std::cout << "Instruction parallelism test:\n---" << std::endl;
	double add1_time, mul1_time, div1_time;
	ArithmeticInstructionParallelismSingleOperation(&add1_time, &mul1_time,
													&div1_time);
	double add2_time, mul2_time, div2_time;
	ArithmeticInstructionParallelismTwoOperations(&add2_time, &mul2_time,
												  &div2_time);
	double add3_time, mul3_time, div3_time;
	ArithmeticInstructionParallelismThreeOperations(&add3_time, &mul3_time,
													&div3_time);
	double add4_time, mul4_time, div4_time;
	ArithmeticInstructionParallelismFourOperations(&add4_time, &mul4_time,
												   &div4_time);
	double add5_time, mul5_time, div5_time;
	ArithmeticInstructionParallelismFiveOperations(&add5_time, &mul5_time,
												   &div5_time);
	double add6_time, mul6_time, div6_time;
	ArithmeticInstructionParallelismSixOperations(&add6_time, &mul6_time,
												  &div6_time);
	double add7_time, mul7_time, div7_time;
	ArithmeticInstructionParallelismSevenOperations(&add7_time, &mul7_time,
													&div7_time);
	double add8_time, mul8_time, div8_time;
	ArithmeticInstructionParallelismEightOperations(&add8_time, &mul8_time,
													&div8_time);
	double add9_time, mul9_time, div9_time;
	ArithmeticInstructionParallelismNineOperations(&add9_time, &mul9_time,
												   &div9_time);
	double add10_time, mul10_time, div10_time;
	ArithmeticInstructionParallelismTenOperations(&add10_time, &mul10_time,
												  &div10_time);

	double add20_time, mul20_time, div20_time;
	ArithmeticInstructionParallelismTwentyOperations(&add20_time, &mul20_time,
													 &div20_time);
	std::cout << "Times with 1 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add1_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul1_time << std::endl;
	std::cout << "Division:\t\t" << div1_time << "\n" << std::endl;

	std::cout << "Times with 2 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add2_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul2_time << std::endl;
	std::cout << "Division:\t\t" << div2_time << "\n" << std::endl;

	std::cout << "Times with 3 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add3_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul3_time << std::endl;
	std::cout << "Division:\t\t" << div3_time << "\n" << std::endl;

	std::cout << "Times with 4 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add4_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul4_time << std::endl;
	std::cout << "Division:\t\t" << div4_time << "\n" << std::endl;

	std::cout << "Times with 5 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add5_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul5_time << std::endl;
	std::cout << "Division:\t\t" << div5_time << "\n" << std::endl;

	std::cout << "Times with 6 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add6_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul6_time << std::endl;
	std::cout << "Division:\t\t" << div6_time << "\n" << std::endl;

	std::cout << "Times with 7 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add7_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul7_time << std::endl;
	std::cout << "Division:\t\t" << div7_time << "\n" << std::endl;

	std::cout << "Times with 8 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add8_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul8_time << std::endl;
	std::cout << "Division:\t\t" << div8_time << "\n" << std::endl;

	std::cout << "Times with 9 operation per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add9_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul9_time << std::endl;
	std::cout << "Division:\t\t" << div9_time << "\n" << std::endl;

	std::cout << "Times with 10 operations per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add10_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul10_time << std::endl;
	std::cout << "Division:\t\t" << div10_time << std::endl;

	std::cout << "Times with 20 operations per math expression: " << std::endl;
	std::cout << "Addition:\t\t" << add20_time << std::endl;
	std::cout << "Multiplication:\t\t" << mul20_time << std::endl;
	std::cout << "Division:\t\t" << div20_time << std::endl;

}

void OperationsBenchmark::ArithmeticInstructionParallelismSingleOperation(double* add_time,
															   double* mul_time, double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a1, a2, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a1, a2, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a1, a2, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}


void OperationsBenchmark::ArithmeticInstructionParallelismTenOperations(double* add_time,
															 double* mul_time, double* div_time) {
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
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front(), * a8 = &a8v.front(), * a9 = &a9v.front(),
			* a10 = &a10v.front(), * a11 = &a11v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
					a8[i] + a9[i] + a10[i] + a11[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i]*
					a8[i]*a9[i]*a10[i]*a11[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i]
					/a8[i]/a9[i]/a10[i]/a11[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}


void
OperationsBenchmark::BenchAll() {
	ArithemticOperations();
	ArithmeticInstructionParallelism();
}

OperationsBenchmark::OperationsBenchmark() {
}


OperationsBenchmark::~OperationsBenchmark() {
}

void OperationsBenchmark::ArithmeticInstructionParallelismTwentyOperations(double* addTime, double* mulTime,
																double* divTime) {
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
			a10v(length), a11v(length), a12v(length), a13v(length), a14v(length),
			a15v(length), a16v(length), a17v(length), a18v(length), a19v(length), a20v(
			length), a21v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front(), * a8 = &a8v.front(), * a9 = &a9v.front(),
			* a10 = &a10v.front(), * a11 = &a11v.front(), * a12 = &a12v.front(), * a13 = &a13v.front(), * a14 = &a14v.front(),
			* a15 = &a15v.front(), * a16 = &a16v.front(), * a17 = &a17v.front(), * a18 = &a18v.front(), * a19 = &a19v.front(),
			* a20 = &a20v.front(), * a21 = &a21v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);
	ResetArrays(length, a12, a13, ignoreit);
	ResetArrays(length, a14, a15, ignoreit);
	ResetArrays(length, a16, a17, ignoreit);
	ResetArrays(length, a18, a19, ignoreit);
	ResetArrays(length, a20, a21, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
					a8[i] + a9[i] + a10[i] + a11[i] + a12[i] + a13[i] + a14[i] + a15[i] + a16[i] +
					a17[i] +
					a18[i] + a19[i] + a20[i] + a21[i];
		}
	}
	sw.Stop();
	*addTime = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);
	ResetArrays(length, a12, a13, ignoreit);
	ResetArrays(length, a14, a15, ignoreit);
	ResetArrays(length, a16, a17, ignoreit);
	ResetArrays(length, a18, a19, ignoreit);
	ResetArrays(length, a20, a21, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i]*
					a8[i]*a9[i]*a10[i]*a11[i]*a12[i]*a13[i]*a14[i]*a15[i]*a16[i]*a17[i]*
					a18[i]*a19[i]*a20[i]*a21[i];
		}
	}
	sw.Stop();
	*mulTime = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArrays(length, a10, a11, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i]
					/a8[i]/a9[i]/a10[i]/a11[i]/a12[i]/a13[i]/a14[i]/a15[i]/a16[i]/a17[i]
					/a18[i]/a19[i]/a20[i]/a21[i];
		}
	}
	sw.Stop();
	*divTime = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismTwoOperations(double* add_time, double* mul_time,
															 double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismThreeOperations(double* add_time, double* mul_time,
															   double* div_time) {
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
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(),
			* a4 = &a4v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArray(length, a4, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArray(length, a4, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArray(length, a4, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismFourOperations(double* add_time, double* mul_time,
															  double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
			a5v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismFiveOperations(double* add_time, double* mul_time,
															  double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
			a5v(length), a6v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArray(length, a6, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArray(length, a6, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArray(length, a6, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismSixOperations(double* add_time, double* mul_time,
															 double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
			a5v(length), a6v(length), a7v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismSevenOperations(double* add_time, double* mul_time,
															   double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
			a5v(length), a6v(length), a7v(length), a8v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front(), * a8 = &a8v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArray(length, a8, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
					a8[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArray(length, a8, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i]*
					a8[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArray(length, a8, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i]
					/a8[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismEightOperations(double* add_time, double* mul_time,
															   double* div_time) {
	const int from = 0;
	const int to = 1;
	std::random_device rd;
	std::mt19937 mt(rd());
	std::uniform_real_distribution<double> dist(from, to);
	StopWatch sw;

	const int length = 512;
	const int loops = 1e5;
	std::vector<double> a0v(length), a1v(length), a2v(length), a3v(length), a4v(length),
			a5v(length), a6v(length), a7v(length), a8v(length), a9v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front(), * a8 = &a8v.front(), * a9 = &a9v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
					a8[i] + a9[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i]*
					a8[i]*a9[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i]
					/a8[i]/a9[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}

void OperationsBenchmark::ArithmeticInstructionParallelismNineOperations(double* add_time, double* mul_time,
															  double* div_time) {
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
			a10v(length);
	double* a0 = &a0v.front(), * a1 = &a1v.front(), * a2 = &a2v.front(), * a3 = &a3v.front(), * a4 = &a4v.front(),
			* a5 = &a5v.front(), * a6 = &a6v.front(), * a7 = &a7v.front(), * a8 = &a8v.front(), * a9 = &a9v.front(),
			* a10 = &a10v.front();
	auto ignoreit = 0.0;

	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArray(length, a10, ignoreit);
	sw.Start();
	for (size_t l = 0; l < loops; l++) {
		// MSVC cannot vectorize this loop (message 1300).
		// However, if this loop will not be nested in, autovectorization will
		//happen. Same condition apply for mul/div/rcp loops

		// ICL does not have this issue, but to provide both vectorized and
		// nonvectorized comparison i specifically disabled vectorization in
		// this method
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i] + a2[i] + a3[i] + a4[i] + a5[i] + a6[i] + a7[i] +
					a8[i] + a9[i] + a10[i];
		}
	}
	sw.Stop();
	*add_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArray(length, a10, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]*a2[i]*a3[i]*a4[i]*a5[i]*a6[i]*a7[i]*
					a8[i]*a9[i]*a10[i];
		}
	}
	sw.Stop();
	*mul_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
	ResetArrays(length, a0, a1, ignoreit);
	ResetArrays(length, a2, a3, ignoreit);
	ResetArrays(length, a4, a5, ignoreit);
	ResetArrays(length, a6, a7, ignoreit);
	ResetArrays(length, a8, a9, ignoreit);
	ResetArray(length, a10, ignoreit);

	sw.Start();
	for (size_t l = 0; l < loops; l++) {
#pragma novector
#pragma loop( no_vector )
		for (int i = 0; i < length; ++i) {
			a0[i] = a1[i]/a2[i]/a3[i]/a4[i]/a5[i]/a6[i]/a7[i]
					/a8[i]/a9[i]/a10[i];
		}
	}
	sw.Stop();
	*div_time = sw.EllapsedTime();
	ignoreit -= a0[((int) (dist(mt))%static_cast<int>(length))];
}