// Copyright 2020 Gysarova Daria
#include <mpi.h>
#include <vector>
#include <random>
#include<algorithm>
#include <ctime>
#include <limits>
#include "../../../modules/task_1/gusarova_d_min_in_matrix/min_in_matr.h"


std::vector<std::vector<int>> RandMatr(int str, int st) {
    std::vector<std::vector<int>> B(str*st);
    std::vector<int> tmp(st);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < str; i++) {
        tmp.clear();
        for (int j = 0; j < st; j++) {
            tmp.push_back(gen());
        }
        B[i] = tmp;
    }
    return B;
}
std::vector<int> MatrInArr(std::vector<std::vector<int>> B, int str, int st) {
    std::vector<int> A(st*str);
    for (int i = 0; i < str; i++) {
        for (int j = 0; j < st; j++) {
            A[i*st+j] = B[i][j];
        }
    }
    return A;
}
int MinInArr(std::vector<int> A) {
    return *std::min_element(A.begin(), A.end());
}
int ParMin(std::vector<int> A) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int ost = A.size() % size;
    int k = A.size() / size;
    int l = 0;
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Send(A.data() + i * k + ost, k, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        l = k + ost;
    } else {
        l = k;
    }
    std::vector<int> LocAr(l);
    if (rank == 0) {
        LocAr = std::vector<int>(A.begin(), A.begin() + k + ost);
    } else {
        MPI_Status status;
        MPI_Recv(&LocAr[0], k, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    int locmin = 0;
    if (LocAr.size() > 0) {
        locmin = *std::min_element(LocAr.begin(), LocAr.end());
    } else {
        locmin = std::numeric_limits<int>::max();
    }
    int min = 0;
    MPI_Reduce(&locmin, &min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    return min;
}
