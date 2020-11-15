// Copyright 2020 Zaikin Ilya

#ifndef MODULES_TASK_1_ZAIKIN_I_MONTECARLO_MONTECARLO_H_
#define MODULES_TASK_1_ZAIKIN_I_MONTECARLO_MONTECARLO_H_

#include <mpi.h>
#include <math.h>
#include <iostream>


float MonteCarloNotParallelMethod(double b, double a, int n, double(*func)(double));
float MonteCarloParallelMethod(double b, double a, int n, double(*func)(double));


#endif  // MODULES_TASK_1_ZAIKIN_I_MONTECARLO_MONTECARLO_H_
