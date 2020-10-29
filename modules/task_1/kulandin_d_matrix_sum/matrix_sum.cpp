// Copyright 2020 Kulandin Denis
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include "../../../modules/task_1/kulandin_d_matrix_sum/matrix_sum.h"

std::vector<int> randomVector(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> ans(sz);
    for (int i = 0; i < sz; ++i) {
        ans[i] = gen();
    }
    return ans;
}

int64_t getParallelSum(const std::vector<int> &a, int n) {
    int rank, size;
    int64_t globalSum = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const int delta = n / size;
    const int rem = n % size;
    if (rank == 0) {
        for (int j = 1; j < size; ++j) {
            MPI_Send(&a[0] + rem + j * delta, delta, MPI_INT, j, 0, MPI_COMM_WORLD);
        }
    }
    std::vector<int> local(rank == 0 ? rem + delta : delta, 0);
    if (rank == 0) {
        local = std::vector<int>(a.begin(), a.begin() + delta + rem);
    } else {
        MPI_Status status;
        MPI_Recv(&local[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    int64_t localSum = getSequentialSum(local);
    MPI_Reduce(&localSum, &globalSum, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return globalSum;
}

int64_t getSequentialSum(const std::vector<int> &vec) {
    return std::accumulate(vec.begin(), vec.end(), 0ll);
}
