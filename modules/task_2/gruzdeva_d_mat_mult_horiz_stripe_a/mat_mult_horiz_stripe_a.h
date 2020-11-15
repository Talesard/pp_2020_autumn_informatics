// Copyright 2020 Gruzdeva Diana
#ifndef MODULES_TASK_2_GRUZDEVA_D_MAT_MULT_HORIZ_STRIPE_A_MAT_MULT_HORIZ_STRIPE_A_H_
#define MODULES_TASK_2_GRUZDEVA_D_MAT_MULT_HORIZ_STRIPE_A_MAT_MULT_HORIZ_STRIPE_A_H_

#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>

std::vector<double> getRandomMatrix(int width, int height, time_t seed);
std::vector<double> getSequentialMultiplication(std::vector<double> matrixA, std::vector<double> matrixB,
                                    int matrixARows, int matrixACols, int matrixBCols);
std::vector<double> getParalleMultiplication(std::vector<double> matrixA, std::vector<double> matrixB,
                                    int matrixARows, int matrixACols, int matrixBRows, int matrixBCols);

#endif  // MODULES_TASK_2_GRUZDEVA_D_MAT_MULT_HORIZ_STRIPE_A_MAT_MULT_HORIZ_STRIPE_A_H_
