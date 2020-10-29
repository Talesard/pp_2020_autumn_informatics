// Copyright 2020 Rustamov Azer
#ifndef MODULES_TASK_1_RUSTAMOV_A_MATRIX_COLS_SUM_MATRIX_COLS_SUM_H_
#define MODULES_TASK_1_RUSTAMOV_A_MATRIX_COLS_SUM_MATRIX_COLS_SUM_H_
#include <vector>

using Matrix = std::vector<int>;

Matrix RandomMatrix(int rows, int cols);
Matrix SequentialColsSum(Matrix matrix, int rows, int cols);
Matrix ParallelColsSum(Matrix matrix, int rows, int cols);

#endif  // MODULES_TASK_1_RUSTAMOV_A_MATRIX_COLS_SUM_MATRIX_COLS_SUM_H_
