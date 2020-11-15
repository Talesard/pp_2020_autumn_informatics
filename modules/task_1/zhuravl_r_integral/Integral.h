// Copyright 2020 Zhuravlev Roman
#ifndef MODULES_TASK_1_ZHURAVL_R_INTEGRAL_INTEGRAL_H_
#define MODULES_TASK_1_ZHURAVL_R_INTEGRAL_INTEGRAL_H_

#include <functional>

double getParallelIntegration(const std::function<double(double)>& f, double a, double b, int n);
double getSequentialIntegration(const std::function<double(double)>& f, double a, double b, int n);

#endif  // MODULES_TASK_1_ZHURAVL_R_INTEGRAL_INTEGRAL_H_
