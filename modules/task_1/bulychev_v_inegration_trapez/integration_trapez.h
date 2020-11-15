// Copyright 2020 Bulychev Vladislav
#ifndef MODULES_TASK_1_BULYCHEV_V_INEGRATION_TRAPEZ_INTEGRATION_TRAPEZ_H_
#define MODULES_TASK_1_BULYCHEV_V_INEGRATION_TRAPEZ_INTEGRATION_TRAPEZ_H_

#include <functional>

double getIntegral(double a, double b, int n, const std::function<double(double)>& f);
double getParallelIntegral(double a, double b, int n, const std::function<double(double)>& f);

#endif  // MODULES_TASK_1_BULYCHEV_V_INEGRATION_TRAPEZ_INTEGRATION_TRAPEZ_H_
