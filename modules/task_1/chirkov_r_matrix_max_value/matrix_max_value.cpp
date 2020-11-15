// Copyright 2020 Chirkov Roman
#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include "../../../modules/task_1/chirkov_r_matrix_max_value/matrix_max_value.h"

std::vector<int> generateMatrix(int rows, int cols) {
    std::vector<int> matrix;
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < rows * cols; i++) {
        matrix.push_back(gen());
    }
    return matrix;
}


int sequentialFind(std::vector<int> matrix, int rows, int cols) {
    int max = matrix[0];
    for (int i = 0; i < rows * cols; i++) {
        if (matrix[i] > max) {
            max = matrix[i];
        }
    }
    return max;
}

int parallelFind(std::vector<int> matrix, int rows, int cols) {
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *sendcounts = new int[size];
    for (int i = 0; i < size; i++) {
        sendcounts[i] = rows * cols / size;
    }
    sendcounts[size-1] += rows * cols % size;

    int *displs = new int[size];
    int sum = 0;
    for (int i = 0; i < size; i++) {
        displs[i] = sum;
        sum+=sendcounts[i];
    }

    int recvcount = rows * cols / size;
    if (rank == size - 1) {
        recvcount+=rows * cols % size;
    }

    int* recvbuf = new int[recvcount];

    MPI_Scatterv(rank == 0 ? &matrix[0] : 0, sendcounts, displs, MPI_INT,
        recvbuf, recvcount, MPI_INT, 0, MPI_COMM_WORLD);

    int partMax = -2147483648;
    for (int i = 0; i < recvcount; i++) {
        if (recvbuf[i] > partMax) {
            partMax = recvbuf[i];
        }
    }

    delete[] sendcounts;
    delete[] displs;
    delete[] recvbuf;

    int globalMax;
    MPI_Allreduce(&partMax, &globalMax, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    return globalMax;
}
