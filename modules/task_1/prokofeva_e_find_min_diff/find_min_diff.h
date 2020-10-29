// Copyright 2020 Prokofeva Elizaveta
#ifndef MODULES_TASK_1_PROKOFEVA_E_FIND_MIN_DIFF_FIND_MIN_DIFF_H_
#define MODULES_TASK_1_PROKOFEVA_E_FIND_MIN_DIFF_FIND_MIN_DIFF_H_

#include <vector>
#include <utility>

using std::vector;
using std::pair;

vector<int> getRandomVector(int size);

pair<int, int> getParallelOperations(const vector<int>& vec);
pair<int, int> getSequentialOperations(const vector<int>& vec);


void pairCompare(int* input, int* output, int* lenght, MPI_Datatype* dtype);

#endif  // MODULES_TASK_1_PROKOFEVA_E_FIND_MIN_DIFF_FIND_MIN_DIFF_H_
