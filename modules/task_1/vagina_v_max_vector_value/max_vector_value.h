// Copyright 2020 Vagina Valeria
#ifndef MODULES_TASK_1_VAGINA_V_MAX_VECTOR_VALUE_MAX_VECTOR_VALUE_H_
#define MODULES_TASK_1_VAGINA_V_MAX_VECTOR_VALUE_MAX_VECTOR_VALUE_H_

#include <mpi.h>
#include <vector>

std::vector<int> getRandomVector(int size);

int maxVectorValue(std::vector<int> vec);
int maxVectorValueParallel(std::vector<int> globalVector);

#endif  // MODULES_TASK_1_VAGINA_V_MAX_VECTOR_VALUE_MAX_VECTOR_VALUE_H_
