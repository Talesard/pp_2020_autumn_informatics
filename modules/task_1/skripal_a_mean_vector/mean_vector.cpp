// Copyright 2020 Skripal Andrey
#include<mpi.h>
#include<vector>
#include <ctime>
#include<random>
#include<iostream>


std::vector<int> createRandomVector(int n) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(n);
    for (int i = 0; i < n; i++) {
        vec[i] = gen()%100;
    }
    return vec;
}

double getSequentialMean(std::vector<int> vec) {
    int size = vec.size();
    double localSum = 0;
    for (int i = 0; i < size; i++) {
        localSum += vec[i];
    }
    double localMean = localSum / size;
    return localMean;
}

double getSequentialSum(std::vector<int> vec) {
    int size = vec.size();
    double localSum = 0;
    for (int i = 0; i < size; i++) {
        localSum += vec[i];
    }
    return localSum;
}

double getMeanVector(std::vector<int> vec, int n) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int delta = n / size;
    int k = n - delta * size;
    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&vec[0] + proc * delta + k, delta, MPI_INT, proc, 0, MPI_COMM_WORLD);
        }
    }
    std::vector<int> localVec;
    if (rank == 0) localVec = std::vector<int>(delta + k);
    else
         localVec = std::vector<int>(delta);

    if (rank == 0) {
         localVec = std::vector<int>(vec.begin(), vec.begin() + delta+k);
    } else {
        MPI_Status status;
        MPI_Recv(&localVec[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    double globalMean = 0;
    double localSum = getSequentialSum(localVec);
    MPI_Reduce(&localSum, &globalMean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    globalMean /= n;
    return globalMean;
}
