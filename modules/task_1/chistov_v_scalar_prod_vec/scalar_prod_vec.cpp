// Copyright 2020 Chistov Vladimir

#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include "../../../modules/task_1/chistov_v_scalar_prod_vec/scalar_prod_vec.h"

std::vector<int> GenRandVec(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(size);
    for (int i = 0; i < size; i++) {
        vec[i] = gen() % 100 + 1;
    }
    return vec;
}

int CalcSumSequential(std::vector<int> a, std::vector<int> b) {
    int size = a.size();
    if (size == 0) {
        return 0;
    }

    int sum = 0;
    for (int i = 0; i < size; i++) {
        sum += a[i] * b[i];
    }
    return sum;
}

int CalcSumParallel(std::vector<int> a, std::vector<int> b) {
    int ProcNum, ProcRank;
    int n = a.size();

    if (n == 0) {
        return 0;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int delta = n / ProcNum;
    int gamma = n % ProcNum;
    std::vector<int> loc_a(delta+gamma);
    std::vector<int> loc_b(delta+gamma);

    if (ProcRank == 0) {
        for (int i = 1; i < ProcNum; i++) {
            MPI_Send(&a[0] + gamma + i * delta, delta, MPI_INT, i, 0, MPI_COMM_WORLD);
            MPI_Send(&b[0] + gamma + i * delta, delta, MPI_INT, i, 1, MPI_COMM_WORLD);
        }
    }

    if (ProcRank == 0) {
        loc_a = std::vector<int>(a.begin(), a.begin() + delta + gamma);
        loc_b = std::vector<int>(b.begin(), b.begin() + delta + gamma);
    } else {
        MPI_Status status;
        MPI_Recv(&loc_a[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&loc_b[0], delta, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
    }

    int sum = 0;
    int loc_sum = 0;
    loc_sum = CalcSumSequential(loc_a, loc_b);

    MPI_Reduce(&loc_sum, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return sum;
}
