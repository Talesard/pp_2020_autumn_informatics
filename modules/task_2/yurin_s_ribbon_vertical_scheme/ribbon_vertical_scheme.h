// Copyright 2020 Yurin Stanislav
#ifndef MODULES_TASK_2_YURIN_S_RIBBON_VERTICAL_SCHEME_RIBBON_VERTICAL_SCHEME_H_
#define MODULES_TASK_2_YURIN_S_RIBBON_VERTICAL_SCHEME_RIBBON_VERTICAL_SCHEME_H_

#include <vector>
#include <string>

unsigned int get_random_time();
std::vector<int> getRandomVector(int  sz);
std::vector<int> getRandomMatrix(int number_of_rows, int number_of_columns);
std::vector<int> getSequentialResultVector(std::vector<int> matrix, int number_of_matrix_rows,
                                            int number_of_matrix_columns, std::vector<int> factor_vector);
std::vector<int> getParallelResultVector(std::vector<int> matrix, int number_of_matrix_rows,
                                            int number_of_matrix_columns, std::vector<int> factor_vector);

#endif  // MODULES_TASK_2_YURIN_S_RIBBON_VERTICAL_SCHEME_RIBBON_VERTICAL_SCHEME_H_
