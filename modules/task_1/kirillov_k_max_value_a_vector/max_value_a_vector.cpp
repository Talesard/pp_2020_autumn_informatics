// Copyright 2020 Kirillov Konstantin
#include <mpi.h>
#include <vector>
#include <ctime>
#include <random>
#include <algorithm>
#include "../../../modules/task_1/kirillov_k_max_value_a_vector/max_value_a_vector.h"

std::vector<int> getRandomVector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int>vec(size);
    for (int i = 0; i < size; i++) {
        vec[i] = gen() % 100;
    }
    return vec;
}
int getSequentialMax(std::vector<int> vec) {
    return *std::max_element(vec.begin(), vec.end());
}
int getParallelMax(std::vector<int> vec, int size_vec) {
    int procNum, procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int delta = size_vec / procNum;
    int tmp = size_vec - delta*procNum;
    if (procRank == 0) {
        for (int proc = 1; proc < procNum; proc++) {
            MPI_Send(&vec[0] + tmp + proc * delta, delta, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }
    int local_size;
    if (procRank == 0) local_size = tmp + delta;
    else
    local_size = delta;
    std::vector<int>local(local_size);
    if (procRank == 0) {
        local = std::vector<int>(vec.begin(), vec.begin() + tmp + delta);
    } else {
        MPI_Status status;
        MPI_Recv(&local[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    int global_Max = 0;
    int local_Max = getSequentialMax(local);
    MPI_Reduce(&local_Max, &global_Max, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    return global_Max;
}
