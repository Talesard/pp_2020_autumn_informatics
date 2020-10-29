// Copyright 2020 Loganov Andrei
#include <mpi.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include <iostream>
#include <algorithm>
#include "../../../../modules/task_1/loganov_a_max_element_matrix/matrix_sum_stolb.h"

std::vector<int> getSequentialMaxStolb(std::vector<int>matrixvec, int cols, int rows) {
    std::vector<int> result(cols);
    for (int i = 0; i < cols; i++) {
        result[i] = (*std::max_element(matrixvec.begin() + i * rows,
            matrixvec.begin() + (i + 1) * rows));
    }
    return result;
}
std::vector<int> randomMatrix(int** Matrix, int cols, int rows) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Matrix[i][j] = gen() % 100;
        }
    }

    std::vector<int> rez(rows * cols);
    for (int j = 0; j < cols; j++) {
        for (int i = 0; i < rows; i++) {
            rez[i + rows * j] = Matrix[i][j];
        }
    }
    return rez;
}
int getSequentialMax(const std::vector<int>& vec) {
    return *std::max_element(vec.begin(), vec.end());
}
std::vector<int> getParallelMaxElements(const std::vector<int>& a, int cols, int rows) {
    int rank, size;
    std::vector<int> glmaxs(cols, 0);
    std::vector<int> glmaxsk(cols, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int delta = cols / size;
    int ost = cols % size;
    std::vector<int> local(rows);
    std::vector<int> deltapoc(delta * rows);
    if (rank == 0) {
        for (int j = 1; j < size; ++j) {
            MPI_Send(&a[0] + ost * rows + j * delta * rows, delta * rows, MPI_INT, j, 0, MPI_COMM_WORLD);
        }
    }
    if (rank == 0) {
        for (int i = 0; i < delta + ost; i++) {
            for (int j = i * rows; j < i * rows + rows; j++) {
                local[j - i * rows] = a[j];
            }
            glmaxs[i] = getSequentialMax(local);
        }
    } else {
        MPI_Status status;
        MPI_Recv(&deltapoc[0], delta * rows, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        for (int i = 0; i < delta; i++) {
            for (int j = rows * i; j < rows * i + rows; j++) {
                local[j - i * rows] = deltapoc[j];
            }
            glmaxsk[i] = getSequentialMax(local);
        }
        MPI_Send(&glmaxsk[0], delta, MPI_INT, 0, 7, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Status status;
            MPI_Recv(&glmaxs[i * delta + ost], delta, MPI_INT, i, 7, MPI_COMM_WORLD, &status);
        }
    }
    return glmaxs;
}
