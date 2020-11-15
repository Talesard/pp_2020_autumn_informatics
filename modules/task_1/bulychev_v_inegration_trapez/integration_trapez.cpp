// Copyright 2020 Bulychev Vladislav
#include <mpi.h>
#include <cmath>
#include <functional>
#include "../../../modules/task_1/bulychev_v_inegration_trapez/integration_trapez.h"

double getIntegral(double a, double b, int n, const std::function<double(double)>& f) {
    if (a > b) {
        throw "Error";
    }
    if (n < 1) {
        throw "Error";
    }
    double result = 0;
    double h;
    h = (b - a) / n;
    if (0 == n) {
        return result;
    }
    for (int i = 1; i < n; ++i) {
        result += f(a + i * h);
    }
    result += (f(a) + f(b)) * 0.5;
    result *= h;
    return result;
}

double getParallelIntegral(double a, double b, int n, const std::function<double(double)>& f) {
    if (a > b) {
        throw "Error";
    }
    if (n < 1) {
        throw "Error";
    }
    double result = 0;
    double h;
    h = (b - a) / n;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (size == 1) {
        return getIntegral(a, b, n, f);
    }
    double l_result = 0;
    double l_n = n / size;
    double l_a = a + h * rank * l_n;
    double l_b = a + rank * h * l_n + h * l_n;
    if (rank != size - 1) {
        l_result = getIntegral(l_a, l_b, l_n, f);
    } else {
        l_result = getIntegral(l_a, b, n - (n / size) * (size - 1.0), f);
    }
    MPI_Reduce(&l_result, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return result;
}
