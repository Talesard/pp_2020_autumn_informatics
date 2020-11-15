// Copyright 2020 Kulandin Denis
#ifndef MODULES_TASK_3_KULANDIN_D_QUICK_SORT_SIMPLE_QUICK_SORT_H_
#define MODULES_TASK_3_KULANDIN_D_QUICK_SORT_SIMPLE_QUICK_SORT_H_
#include <vector>

std::vector<int> randomVector(int sz);
void division(std::vector<int>* a, int l, int r, int* t);
void quickSort(std::vector<int>* a, int l, int r);
std::vector<int> parallelQuickSort(const std::vector<int>& a);

#endif  // MODULES_TASK_3_KULANDIN_D_QUICK_SORT_SIMPLE_QUICK_SORT_H_
