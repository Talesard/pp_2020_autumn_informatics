// Copyright 2020 Galkin Vladislav
#include <mpi.h>
#include <vector>
#include<algorithm>
#include <random>
#include "../../../../modules/task_1/galkin_v_min_val/min_values_of_strings.h"

std::vector<std::vector<int>> createRandomMatrix(int n, int m) {
    std::vector<std::vector<int>> matrix;
    std::random_device dev;
    std::mt19937 gen(dev());
    for (int i = 0; i < n; i++) {
        std::vector<int> str;
        for (int j = 0; j < m; j++) {
            str.push_back(static_cast<int>(gen() % 1000));
        }
        matrix.push_back(str);
    }
    return matrix;
}

std::vector<int> getMinSequential(std::vector<std::vector<int>> matrix) {
    if (matrix.size() == 0 || matrix[0].size() == 0) throw "Invalid size matrix";
    int n = matrix.size(), m = matrix[0].size();
    std::vector<int> res;
    for (int i = 0; i < n; i++) {
        int min = matrix[i][0];
        for (int j = 0; j < m; j++) {
            if (min > matrix[i][j]) min = matrix[i][j];
        }
        res.push_back(min);
    }
    return res;
}

std::vector<int> getMinParallel(std::vector<std::vector<int>> matrix) {
    if (matrix.size() == 0 || matrix[0].size() == 0) throw "Invalid size matrix";
    int n = matrix.size(), m = matrix[0].size();
    std::vector<int> minval(n);
    std::vector<int> tmp;
    int procnum, rank;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &procnum);
    int t = n / procnum;
    int rem = 0;
    if (t != 0) rem = n % procnum;
    if (rank == 0) {
        if (n > 1) {
            for (int i = 1; i < std::min(n, procnum); i++) {
                for (int j = 0; j < n; j++) {
                    MPI_Send(matrix[j].data(), m, MPI_INT, i, 20 + j, MPI_COMM_WORLD);
                }
            }
        }
        if (t == 0) t++;
        for (int i = 0; i < t; i++) {
            int min = matrix[i][0];
            for (int j = 0; j < m; j++) {
            if (min > matrix[i][j]) min = matrix[i][j];
            }
            minval[i] = min;
        }
        if (n / procnum == 0) t--;
        if (t != 0 && rem != 0) {
            int min = matrix[n - rem][0];
            for (int i = 0; i < m; i++) {
                if (min > matrix[n - rem][i]) min = matrix[n - rem][i];
            }
            minval[n - rem] = min;
        }
        for (int i = 1; i < std::min(n, procnum); i++) {
            if (t != 0) {
                MPI_Recv(&minval[n / procnum * i], n / procnum, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            } else {
                MPI_Recv(&minval[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
            }
        }
        for (int i = 1; i < rem; i++) {
            MPI_Recv(&minval[n - rem + i], 1, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
        }
    } else {
        if (rank < n) {
            for (int i = 0; i < n; i++) {
                MPI_Recv(matrix[i].data(), m, MPI_INT, 0, 20 + i, MPI_COMM_WORLD, &status);
            }
            if (t == 0) t++;
            for (int i = t * rank; i < t * (rank + 1); i++) {
                if (i < n) {
                    int min = matrix[i][0];
                    for (int j = 0; j < m; j++) {
                        if (min > matrix[i][j]) min = matrix[i][j];
                    }
                    tmp.push_back(min);
                }
            }
            if (n / procnum == 0) t--;
            if (t != 0) {
                MPI_Send(tmp.data(), n / procnum, MPI_INT, 0, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(tmp.data(), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
            }
            if (rem > 1 && t != 0) {
                int k = n - rem + rank;
                if (k < n) {
                    int min = matrix[k][0];
                    for (int i = 0; i < m; i++) {
                        if (min > matrix[k][i]) min = matrix[k][i];
                    }
                    MPI_Send(&min, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
                }
            }
        }
    }
    return minval;
}
