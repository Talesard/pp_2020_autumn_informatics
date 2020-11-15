// Copyright 2020 Gruzdeva Diana

#include "../../../modules/task_2/gruzdeva_d_mat_mult_horiz_stripe_a/mat_mult_horiz_stripe_a.h"
#include <random>
#include <vector>
#include <ctime>

std::vector<double> getRandomMatrix(int rows, int cols, time_t seed) {
    std::vector<double> mat(rows * cols);
    std::mt19937 gen(seed);
    std::uniform_real_distribution<double> dis(-10, 10);
    for (int i = 0; i < rows * cols; i++) {
        mat[i] = dis(gen);
    }
    return mat;
}

std::vector<double> getSequentialMultiplication(std::vector<double> matrixA, std::vector<double> matrixB,
                                    int matrixARows, int matrixACols, int matrixBCols) {
    std::vector<double> matrixC(matrixARows * matrixBCols);
    for (int i = 0; i < matrixARows; i++) {
        for (int j = 0; j < matrixBCols; j++) {
            matrixC[i * matrixBCols + j] = 0;
            for (int k = 0; k < matrixACols; k++) {
                matrixC[i * matrixBCols + j] += matrixA[i * matrixACols + k] * matrixB[k * matrixBCols + j];
            }
        }
    }
    return matrixC;
}

std::vector<double> getParalleMultiplication(std::vector<double> matrixA, std::vector<double> matrixB,
                                    int matrixARows, int matrixACols, int matrixBRows, int matrixBCols) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<double> matrixC;
    std::vector<int> stripeSize(size);
    std::vector<int> displs(size);
    int stripeHeight = matrixARows / size;
    for (int i = 0; i < size; i++) {
        stripeSize[i] = stripeHeight * matrixACols;
        displs[i] = i * stripeHeight * matrixACols;
    }
    int remaining = matrixARows - (size - 1) * stripeHeight;
    stripeSize[size - 1] = remaining * matrixACols;
    if (rank == size - 1) {
        stripeHeight = remaining;
    }
    std::vector<double> matrixAStripe(stripeHeight * matrixACols);

    MPI_Scatterv(matrixA.data(), stripeSize.data(), displs.data(), MPI_DOUBLE,
                 matrixAStripe.data(), stripeSize[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    std::vector<double> matrixCStripe(stripeHeight * matrixBCols);
    for (int i = 0; i < stripeHeight; i++) {
        for (int j = 0; j < matrixBCols; j++) {
            double s = 0;
            for (int k = 0; k < matrixACols; k++) {
                s += matrixAStripe[i * matrixACols + k] * matrixB[k * matrixBCols + j];
            }
            matrixCStripe[i * matrixBCols + j] = s;
        }
    }

    if (rank == 0) {
        matrixC.resize(matrixARows * matrixBCols);
    }

    stripeHeight = matrixARows / size;

    for (int i = 0; i < size; i++) {
        stripeSize[i] = stripeHeight * matrixBCols;
        displs[i] = i * stripeHeight * matrixBCols;
    }
    stripeSize[size - 1] = remaining * matrixBCols;

    MPI_Gatherv(matrixCStripe.data(), stripeSize[rank], MPI_DOUBLE, matrixC.data(), stripeSize.data(), displs.data(),
                                                            MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return matrixC;
}
