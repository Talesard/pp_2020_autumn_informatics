// Copyright 2020 Lyutova Tanya
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include"../../modules/task_1/lyutova_t_matrix_max_in_line/matrix_max_in_line.h"

std::vector<int> getRandomMatrix(int rows, int columns) {
    if (rows < 0 || columns < 0)
        throw "Matrix is incorrect";

    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> matrix(rows * columns);
    for (int i = 0; i < rows * columns; i++) {
        matrix[i] = gen() % 100;
    }
    return matrix;
}

std::vector<int> getMaxInLineSequental(const std::vector<int>matrix, int rows, int columns) {
    if (rows <= 0 || columns <= 0)
        throw "Matrix is incorrect";

    std::vector<int> res(rows);
    for (int i = 0; i < rows; i++) {
        res[i] = (*std::max_element(matrix.begin() +  i * columns,
            matrix.begin() + (i + 1) * columns));
    }
    return res;
}

std::vector<int> getMaxInLineParallel(const std::vector<int>matrix, int rows, int columns) {
    if (rows <= 0 || columns <= 0)
        throw "Matrix is incorrect";

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = rows / size;
    const int rem = rows % size;
    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(matrix.data() + rem * columns + proc * delta * columns,
                delta * columns, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }

    int local_size = 0;
    if (rank == 0)
        local_size = (delta + rem) * columns;
    else
        local_size = delta * columns;
    std::vector<int> local_vec(local_size);

    if (rank == 0) {
        local_vec = std::vector<int>(matrix.begin(), matrix.begin() + local_size);
    } else {
        MPI_Status status;
        MPI_Recv(local_vec.data(), local_size,
            MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    int local_rows = delta;
    if (rank == 0)
        local_rows += rem;
    std::vector<int> local_max(local_rows);
    if (local_rows > 0 && columns > 0)
        local_max = getMaxInLineSequental(local_vec, local_rows, columns);

    if (rank != 0) {
        MPI_Send(local_max.data(), delta, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    std::vector<int> res(rows);

    if (rank == 0) {
        for (int i = 0; i < rem + delta; i++)
            res[i] = local_max[i];
        for (int i = 1; i < size; i++) {
            MPI_Status status;
            MPI_Recv(res.data() + rem + i * delta, delta,
                MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }
    }
    return res;
}
