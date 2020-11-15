// Copyright 2020 Belik Julia
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <ctime>
#include <vector>
#include <random>
#include <algorithm>
#include "../../../modules/task_2/belik_j_ribbon_horiz/RibbonHoriz.h"
std::vector<double> Get_sequestional(std::vector<double> matrix, int rows, int cols, std::vector<double> vect) {
    std::vector<double> vec(cols);
    for (int i = 0; i < cols; i++) {
        vec[i] = 0;
    }
    for (int i = 0; i < cols; i++) {
        vec[i] = skal(i, matrix, rows, vect);
    }
    return vec;
}
double skal(int i, std::vector<double> matrix, int rows, std::vector<double> vect) {
    double sum = 0;
    for (int j = 0; j < rows; j++)
        sum += matrix[j + i * rows] * vect[j];
    return sum;
}
std::vector<double> Get_parallel(std::vector<double> matrix, int rows, int cols, std::vector<double> vect) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int countint = cols / size;
    int rem = cols % size;
    std::vector<double> locmatrix(countint * rows);
    MPI_Status Status;
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Send(matrix.data() + i * countint * rows, countint * rows, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        for (int i = 0; i < countint * rows; i++)
            locmatrix[i] = matrix[i];
    } else {
        MPI_Recv(locmatrix.data(), countint * rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
    }
    std::vector<double> mult(rank == 0 ? cols : countint);
    for (int i = 0; i < countint; i++) {
        mult[i] = skal(i, locmatrix, rows, vect);
    }
    std::vector<double> rmult(cols);
    if (rank == 0) {
        for (int i = 0; i < countint; i++)
            rmult[i] = mult[i];
        for (int i = 1; i < size; i++) {
            MPI_Recv(rmult.data() + i * countint, countint, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &Status);
        }
    } else {
        MPI_Send(mult.data(), countint, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    if (rank == 0) {
        if (rem != 0) {
            for (int j = cols - rem; j < cols; j++)
                rmult[j] = skal(j, matrix, rows, vect);
        }
    }
    return rmult;
}
std::vector<double> Matrix(int rows, int cols) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<double> arr(cols * rows);
    for (int i = 0; i < cols * rows; i++)
        arr[i] = gen() % 100;
    return arr;
}
std::vector<double> Vector(int n) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<double> vec(n);
    for (int i = 0; i < n; i++)
        vec[i] = gen() % 100;
    return vec;
}
