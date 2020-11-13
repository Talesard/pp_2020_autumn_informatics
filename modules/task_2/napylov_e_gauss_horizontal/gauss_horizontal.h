// Copyright 2020 Napylov Evgenii
#ifndef MODULES_TASK_2_NAPYLOV_E_GAUSS_HORIZONTAL_GAUSS_HORIZONTAL_H_
#define MODULES_TASK_2_NAPYLOV_E_GAUSS_HORIZONTAL_GAUSS_HORIZONTAL_H_
#include <vector>

void print_vec(std::vector<double> vec);

void print_vec(std::vector<int> vec);

void print_matrix(std::vector<double> vec, int rows, int cols);

std::vector<double> SystemForPerformanceTest(int rowsm, int cols);

std::vector<double> SolveGaussSeq(std::vector<double> sys, int rows, int cols);

bool CheckSolution(std::vector<double> sys, int rows, int cols, std::vector<double> answer, double epsilon);

std::vector<double> SolveGaussParallel(std::vector<double> sys, int rows, int cols);

#endif  // MODULES_TASK_2_NAPYLOV_E_GAUSS_HORIZONTAL_GAUSS_HORIZONTAL_H_
