// Copyright 2020 Belik Julia
#ifndef MODULES_TASK_1_BELIK_J_INTEGRATION_TRAPEZ_INTEGRATIONTRAPEZ_H_
#define MODULES_TASK_1_BELIK_J_INTEGRATION_TRAPEZ_INTEGRATIONTRAPEZ_H_

#include <functional>

double getParallelIntegration(const std::function<double(double)>& f, double a, double b, int n);
double getSequentialIntegration(const std::function<double(double)>& f, double a, double b, int n);

#endif  // MODULES_TASK_1_BELIK_J_INTEGRATION_TRAPEZ_INTEGRATIONTRAPEZ_H_
