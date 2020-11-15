// Copyright 2020 Tareev Daniil
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/tareev_d_col_min/col_min.h"

std::vector<int> randomMatrix(int cols, int rows) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> rmat(cols * rows);
    for (int i = 0; i < cols * rows; i++) {
        rmat[i] = gen() % 100;
    }
    return rmat;
}

std::vector<int> transposeMatrix(std::vector<int> mat, int cols, int rows) {
    std::vector<int> tmat(cols * rows);
    for (int i = 0; i < cols; i++)
        for (int j = 0; j < rows; j++)
            tmat[j + i * rows] = mat[i + j * cols];
    return tmat;
}

std::vector<int> getSequentialOperations(std::vector<int> mat, int cols, int rows) {
    std::vector<int> min_elems(rows);
    for (int i = 0; i < rows; i++) {
        int min_elem = mat[i * cols];
        for (int j = 1; j < cols; j++)
            min_elem = std::min(min_elem, mat[i * cols + j]);
        min_elems[i] = min_elem;
    }
    return min_elems;
}

std::vector<int> getParallelOperations(std::vector<int> global_mat, int cols, int rows) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = rows / size;
    const int epsilon = rows % size;
    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(global_mat.data() + proc * delta * cols + epsilon * cols,
                     delta * cols, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }
    std::vector<int> local_mat(rank == 0 ? (delta + epsilon) * cols : delta * cols);
    if (rank == 0) {
        local_mat = {global_mat.begin(), global_mat.begin() + (delta + epsilon) * cols};
    } else {
        MPI_Status status;
        MPI_Recv(local_mat.data(), delta * cols, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    std::vector<int> result(rows);
    std::vector<int> local_result = getSequentialOperations(local_mat, cols, rank == 0 ? delta + epsilon : delta);
    if (rank == 0) {
        for (int i = 0; i < delta + epsilon; i++)
            result[i] = local_result[i];
        for (int proc = 1; proc < size; proc++) {
            MPI_Status status;
            MPI_Recv(result.data() + epsilon + delta * proc, delta, MPI_INT, proc, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        MPI_Send(local_result.data(), delta, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    return result;
}
