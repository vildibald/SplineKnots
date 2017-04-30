#pragma once

#include <forward_list>

class OperationsBenchmark {
    void
    ResetArray(const int length, double *a, double &ignoreit);

    void
    ResetArrays(const int length, double* a, double* b, double& ignoreit);
public:

    void
    ArithemticOperations();

    void
    ArithmeticInstructionParallelism();

    void
    BenchAll();

    OperationsBenchmark();

    ~OperationsBenchmark();

private:
    void
    ArithmeticInstructionParallelismSingleOperation(double *add_time, double *mul_time,
                                                    double *div_time);

    void
    ArithmeticInstructionParallelismTwoOperations(double *add_time, double *mul_time,
                                                  double *div_time);

    void
    ArithmeticInstructionParallelismThreeOperations(double *add_time, double *mul_time,
                                                    double *div_time);

    void
    ArithmeticInstructionParallelismFourOperations(double *add_time, double *mul_time,
                                                   double *div_time);

    void
    ArithmeticInstructionParallelismFiveOperations(double *add_time, double *mul_time,
                                                   double *div_time);

    void
    ArithmeticInstructionParallelismSixOperations(double *add_time, double *mul_time,
                                                  double *div_time);

    void
    ArithmeticInstructionParallelismSevenOperations(double *add_time, double *mul_time,
                                                    double *div_time);

    void
    ArithmeticInstructionParallelismEightOperations(double *add_time, double *mul_time,
                                                    double *div_time);

    void
    ArithmeticInstructionParallelismNineOperations(double *add_time, double *mul_time,
                                                   double *div_time);

    void
    ArithmeticInstructionParallelismTenOperations(double *add_time, double *mul_time,
                                                  double *div_time);

    void
    ArithmeticInstructionParallelismTwentyOperations(double *addTime, double *mulTime,
                                                     double *divTime);
};