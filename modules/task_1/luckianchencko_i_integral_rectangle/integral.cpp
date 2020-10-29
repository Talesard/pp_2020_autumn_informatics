// Copyright 2020 Luckyanchencko Ivan
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include"../../../modules/task_1/luckianchencko_i_integral_rectangle/integral.h"

double func1(double x) { return 1 / (log(pow(x, 2))+x); }
double func2(double x) { return exp(x) / (1+pow(x, 3)); }
double func3(double x) { return log(x) / (exp(pow(x, 3)) +x); }

double get_Paral_Integral(double(*func)(double), int a, int b) {
    int n = 10, ProcNum, ProcRank;
    double ans = 0.0, res = 0.0;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double h = static_cast<double>(b - a) / static_cast<double>(n);
    for (int i = ProcRank; i < n; i += ProcNum) {
        double x = static_cast<double>(a) + h * static_cast<double>(i) + h / 2;
        res += func(x);
    }
    MPI_Reduce(&res, &ans, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    ans *= h;
    return ans;
}

double get_Integral(double(*func)(double), int a, int b) {
    int n = 10;
    double  h,  res = 0.0;
    h = static_cast<double>(b - a) / static_cast<double>(n);
    for (int i =0; i < n; i++) {
        double x = static_cast<double>(a) + h * static_cast<double>(i) + h/2;
        res += func(x);
    }
    res = res * h;
    return res;
}
