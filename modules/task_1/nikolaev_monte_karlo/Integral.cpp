// Copyright 2020 Nikolaev Denis
#include <mpi.h>
#include <random>
#include <ctime>
#include <numeric>
#include <iostream>
#include <algorithm>
#include <cstring>
#include "../../../modules/task_1/nikolaev_monte_karlo/Integral.h"
#define SIZE 20000

void getOperations(double(*func)(double), double* x, double* y, double* top, double l, double r, double* Spar) {
    const double step = 0.001;
    double maxY = func(l);
    for (double x = l; x <= r; x += step) {
        const double max = func(x);
        if (max > maxY) {
        maxY = max;
        }
    }
    *top = (maxY + 1);
    srand(time(NULL));
    std::mt19937 gen(time(0));
    std::uniform_real_distribution<> urx(l, r);
    std::uniform_real_distribution<> ury1(0, *top);
    bool alreadyThere;
    for (int i = 0; i < SIZE;) {
        alreadyThere = false;
        double newRandomValueX = urx(gen);
        double newRandomValueY = ury1(gen);
        for (int j = 0; j < i; j++) {
            if ((x[j] == newRandomValueX) || (y[j] == newRandomValueY)) {
                alreadyThere = true;
                break;
            }
        }
        if (!alreadyThere) {
            x[i] = newRandomValueX;
            y[i] = newRandomValueY;
            i++;
        }
    }
    *Spar = (r - l) * (*top);
}
double f1(const double x) {
    return (sin(x) / 2);
}
double f2(const double x) {
    return ((2 * pow(x, 2)) / log(3));
}
double f3(const double x) {
    return ((x * exp(1))/5);
}
double calculating_the_integral(double(*func)(double), double* x, double* y, double* Spar) {
    int count = 0;
    for (int i = 0; i < SIZE; i++) {
        if (y[i] <= func(x[i])) {
            count++;
        }
    }
    double integral = ((count * (*Spar)) / static_cast<double>(SIZE));
    return integral;
}
double calculating_the_paral_integral(double(*func)(double), double* x, double* y, double* Spar) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int delta = SIZE / size;
    int ost = SIZE % size;
    int global_res = 0;
    double* local_x;
    local_x = new double[delta];
    double* local_y;
    local_y = new double[delta];
    int global_res_zero = 0;
    if (rank == 0) {
        for (int j = 1; j < size; j++) {
            MPI_Send(&x[0] + ost + j * delta, delta, MPI_DOUBLE, j, 0, MPI_COMM_WORLD);
            MPI_Send(&y[0] + ost + j * delta, delta, MPI_DOUBLE, j, 1, MPI_COMM_WORLD);
        }
    }
    if (rank == 0) {
        for (int i = 0; i < delta + ost; i++) {
            if (y[i] < func(x[i])) {
                global_res_zero++;
            }
        }
    } else {
        MPI_Status status1;
        MPI_Status status2;
        MPI_Recv(&local_x[0], delta, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status1);
        MPI_Recv(&local_y[0], delta, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status2);
    }
    if (delta + ost != SIZE) {
        int local_res = 0;
        for (int i = 0; i < delta; i++) {
            if (local_y[i] < func(local_x[i])) {
                local_res++;
            }
        }
        MPI_Reduce(&local_res, &global_res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    global_res += global_res_zero;
    double integral = ((global_res * (*Spar)) / static_cast<double>(SIZE));
    return integral;
}
