// Copyright 2020 Yurin Stanislav
#ifndef MODULES_TASK_1_YURIN_S_SUM_OF_VALUES_BY_MATRIX_ROWS_SUM_OF_VALUES_BY_MATRIX_ROWS_H_
#define MODULES_TASK_1_YURIN_S_SUM_OF_VALUES_BY_MATRIX_ROWS_SUM_OF_VALUES_BY_MATRIX_ROWS_H_

#include <vector>
#include <string>

std::vector<int> getSequentialSumOfMatrixRows(std::vector<int> matrix,
                                                int number_of_rows);
std::vector<int> getRandomVector(int  sz);
int getSumOfOneRow(std::vector<int> row);
std::vector<int> getSumOfOneProc(std::vector<int> row,
                                    int number_of_elements_in_one_row);
std::vector<int> getParallelSumOfMatrixRows(std::vector<int> matrix,
                                            int number_of_rows);

#endif  // MODULES_TASK_1_YURIN_S_SUM_OF_VALUES_BY_MATRIX_ROWS_SUM_OF_VALUES_BY_MATRIX_ROWS_H_
