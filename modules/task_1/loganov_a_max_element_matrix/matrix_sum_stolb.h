// Copyright 2020 Loganov Andrei
#ifndef  MODULES_TASK_1_LOGANOV_A_MAX_ELEMENT_MATRIX_MATRIX_SUM_STOLB_H_
#define  MODULES_TASK_1_LOGANOV_A_MAX_ELEMENT_MATRIX_MATRIX_SUM_STOLB_H_

#include <vector>
std::vector<int> randomMatrix(int** Matrix, int cols, int rows);
int getSequentialMax(const std::vector<int>& vec);
std::vector<int> getParallelMaxElements(const std::vector<int>& a, int cols, int rows);
std::vector<int> getSequentialMaxStolb(std::vector<int>matrixvec, int cols, int rows);
#endif  //  MODULES_TASK_1_LOGANOV_A_MAX_ELEMENT_MATRIX_MATRIX_SUM_STOLB_H_

