// Copyright 2020 Nikolaev Denis
#pragma once
#ifndef MODULES_TASK_1_NIKOLAEV_MONTE_KARLO_INTEGRAL_H_
#define MODULES_TASK_1_NIKOLAEV_MONTE_KARLO_INTEGRAL_H_

void getOperations(double(*func)(double), double* x, double* y, double* top, double l, double r, double* Spar);
double calculating_the_paral_integral(double(*func)(double), double* x, double* y, double* Spar);
double calculating_the_integral(double(*func)(double), double* x, double* y, double* Spar);
double f1(double x);
double f2(double x);
double f3(double x);

#endif  // MODULES_TASK_1_NIKOLAEV_MONTE_KARLO_INTEGRAL_H_
