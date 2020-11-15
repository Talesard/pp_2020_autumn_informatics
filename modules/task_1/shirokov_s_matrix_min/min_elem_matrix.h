// Copyright 2020 Shirokov Sergey
#ifndef MODULES_TASK_1_SHIROKOV_S_MATRIX_MIN_MIN_ELEM_MATRIX_H_
#define MODULES_TASK_1_SHIROKOV_S_MATRIX_MIN_MIN_ELEM_MATRIX_H_
#include <random>
#include <iostream>
#include <vector>

void matrix_enter(int sizel_f, int sizec_f, int* a[]);
void enter_massive(int sizel_f, int* a);
void recurv_massive(int sizel_f, int* a);
void matrix_print(int sizel_f, int sizec_f, int* a[]);
void massive_ent(int sizel_f, int sizec_f, int* a[], int a1[]);
int min_el(int sizec_f, int a[]);
int maxElementMatrix(int sizel_f, int sizec_f);

#endif  // MODULES_TASK_1_SHIROKOV_S_MATRIX_MIN_MIN_ELEM_MATRIX_H_
