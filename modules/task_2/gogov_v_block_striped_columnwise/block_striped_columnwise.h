// Copyright 2020 Gogov Vladislav
#ifndef MODULES_TASK_2_GOGOV_V_BLOCK_STRIPED_COLUMNWISE_BLOCK_STRIPED_COLUMNWISE_H_
#define MODULES_TASK_2_GOGOV_V_BLOCK_STRIPED_COLUMNWISE_BLOCK_STRIPED_COLUMNWISE_H_
#include <vector>

bool assertVectors(const std::vector<double>& vec1, const std::vector<double>& vec2);

using Matrix = std::vector<double>;

void printVector(const std::vector<double>& vec, int vecSize);
void printMatrix(const std::vector<double>& vec, int rows, int cols);
Matrix createRandomMatrix(int rows, int cols);

void arraysDistribution(int* sendNum, int* sendInd, int n, int procNum, int remain, int cols);
std::vector<double> sumMatrixByCols(const std::vector<double>& vec, int rows, int cols);

std::vector<double> multMatrixByVectorSequential(const Matrix& matirx, int rows, int cols,
                                            const std::vector<double>& vec, int vecSize);
std::vector<double> multMatrixByVectorParallel(const Matrix& matirx, int rows, int cols,
                                            const std::vector<double>& vec, int vecSize);

#endif  // MODULES_TASK_2_GOGOV_V_BLOCK_STRIPED_COLUMNWISE_BLOCK_STRIPED_COLUMNWISE_H_
