// Copyright 2020 Belik Julia
#ifndef MODULES_TASK_2_BELIK_J_RIBBON_HORIZ_RIBBONHORIZ_H_
#define MODULES_TASK_2_BELIK_J_RIBBON_HORIZ_RIBBONHORIZ_H_

#include <vector>

std::vector<double> Get_parallel(std::vector<double> matrix, int rows, int cols, std::vector<double> vector);
std::vector<double> Get_sequestional(std::vector<double> matrix, int rows, int cols, std::vector<double> vector);
std::vector<double> Matrix(int rows, int cols);
std::vector<double> Vector(int n);
double skal(int i, std::vector<double> matrix, int rows, std::vector<double> vector);

#endif  // MODULES_TASK_2_BELIK_J_RIBBON_HORIZ_RIBBONHORIZ_H_
