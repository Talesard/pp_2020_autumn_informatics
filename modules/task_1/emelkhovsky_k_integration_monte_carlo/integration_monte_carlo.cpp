// Copyright 2020 Emelkhovsky Ekaterina
#include <mpi.h>
#include <stdlib.h>
#include <vector>
#include <random>
#include <ctime>
#include <cmath>
#include <exception>
#include "../../../modules/task_1/emelkhovsky_k_integration_monte_carlo/integration_monte_carlo.h"



int funcF(int x) {
    return x + 1;
}


double EasyFunc(int N, int a, int b) {
    double y = 0;
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 1; i < N; i++) {
        int dot = a + gen() % (b - a);
        y += funcF(dot);
    }
    double dop = static_cast<double>(b - a) / N;
    return (dop * y);
}


double ParallelFunc(int N, int a, int b) {
    if (b < a)
        throw "[Oops...] b must be bigger than a!";
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&a, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int count_dot;
    count_dot = N / size;
    if (rank == 0) {
        count_dot += N % size;
    }

    int dot, y = 0;
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 1; i < count_dot; i++) {
        dot = a + gen() % (b - a);
        y += funcF(dot);
    }
    int global_sum;
    MPI_Reduce(&y, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum;
}
