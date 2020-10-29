// Copyright 2020 Rustamov Azer
#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include <numeric>
#include "../../../modules/task_1/rustamov_a_matrix_cols_sum/matrix_cols_sum.h"

Matrix RandomMatrix(int rows, int cols) {
    std::random_device rd;
    std::mt19937 mersenne(rd());
    Matrix matrix(rows * cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = static_cast<int>(mersenne() % 100u);
        }
    }
    return matrix;
}

Matrix SequentialColsSum(Matrix matrix, int rows, int cols) {
    if (rows * cols == 0) {
        throw "Invalid matrix size.";
    }
    Matrix sum(cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            sum[j] += matrix[i * cols + j];
        }
    }
    return sum;
}

Matrix ParallelColsSum(Matrix matrix, int rows, int cols) {
    if (rows * cols == 0) {
        throw "Invalid matrix size";
    }
    if (rows == 1) {
        return matrix;
    }
    if (cols == 1) {
        Matrix sum(1);
        sum[0] = std::accumulate(matrix.begin(), matrix.end(), 0);
        return sum;
    }
    int procNum, procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    const int n = rows / procNum;
    int delta = rows % procNum;
    Matrix global_sum(cols);
    Matrix local_sum(cols);
    if (procRank == 0) {
        for (int i = 1; i < procNum; i++) {
            MPI_Send(&matrix[0] + (delta * cols) + (i * n * cols), n * cols, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
    }
    Matrix local_matrix((n + (procRank == 0 ? delta : 0)) * cols);
    if (procRank == 0) {
        local_matrix = Matrix(matrix.begin(), matrix.begin() + (n + delta) * cols);
    } else {
        MPI_Status status;
        MPI_Recv(&local_matrix[0], n * cols, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    local_sum = SequentialColsSum(local_matrix, n + (procRank == 0 ? delta : 0), cols);
    MPI_Reduce(&local_sum[0], &global_sum[0], cols, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum;
}
