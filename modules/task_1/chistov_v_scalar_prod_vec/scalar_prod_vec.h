// Copyright 2020 Chistov Vladimir

#ifndef MODULES_TASK_1_CHISTOV_V_SCALAR_PROD_VEC_SCALAR_PROD_VEC_H_
#define MODULES_TASK_1_CHISTOV_V_SCALAR_PROD_VEC_SCALAR_PROD_VEC_H_

#include <vector>

std::vector<int> GenRandVec(int  size);
int CalcSumParallel(std::vector<int> a, std::vector<int> b);
int CalcSumSequential(std::vector<int> a, std::vector<int> b);

#endif  // MODULES_TASK_1_CHISTOV_V_SCALAR_PROD_VEC_SCALAR_PROD_VEC_H_
