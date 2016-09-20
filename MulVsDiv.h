#pragma once
#include <forward_list>

class MulVsDiv
{
	void
	ResetArrays(const int length, double* a, double* b, double& ignoreit);

	void
	ResetMatrix(const int rows, const int columnss, double** matrix, 
		double& ignoreit);

	void
	ResetList(const int length, std::forward_list<double>& list, 
		double& ignoreit);
public:
	void
	Loop();

	void
	ArrayAndNumberLoop();

	void
	LoopVectorized();

	void
	DynamicArrayLoop();

	void
	DynamicArrayLoopVectorized();

	void
	DynamicListLoop();

	void
	DynamicArrayAndNumberLoop();

	void
	CsabaDynamicArrayLoop();

	void
	BackwardDependendDynamicArrayLoop();

	void
	ForwardDependendDynamicArrayLoop();
	
	void
	DependendDynamicArrayLoop();
	
	void
	MemCpy();

	void
	ArithmeticInstructionParallelism();
	
	void
	BenchAll();

	MulVsDiv();

	~MulVsDiv();

private:
	void
	ArithmeticInstructionParallelismSingleOperand(double* add_time, double* mul_time, double* div_time);
	void
	ArithmeticInstructionParallelismTenOperands(double* add_time, double* mul_time, double* div_time);
};
