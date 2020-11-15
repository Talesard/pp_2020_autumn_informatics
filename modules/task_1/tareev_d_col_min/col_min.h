// Copyright 2020 Tareev Daniil
#ifndef MODULES_TASK_1_TAREEV_D_COL_MIN_COL_MIN_H_
#define MODULES_TASK_1_TAREEV_D_COL_MIN_COL_MIN_H_

#include <vector>

std::vector<int> randomMatrix(int cols, int rows);
std::vector<int> transposeMatrix(std::vector<int> mat, int cols, int rows);
std::vector<int> getParallelOperations(std::vector<int> global_mat, int cols, int rows);
std::vector<int> getSequentialOperations(std::vector<int> mat, int cols, int rows);

#endif  // MODULES_TASK_1_TAREEV_D_COL_MIN_COL_MIN_H_
