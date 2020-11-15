// Copyright 2020 Molotkova Svetlana
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <cmath>
#include <ctime>
#include <vector>
#include <random>
#include <iomanip>
#include <map>
#include <algorithm>
#include "../../../modules/task_1/molotkova_s_max_in_cols/Max.h"

double findmax(int stolbnum, double** arr, int N) {
    double max = 0;
    max = arr[0][stolbnum];
    for (int i = 0; i < N; i++)
        if (arr[i][stolbnum] > max)
            max = arr[i][stolbnum];
    return max;
}
int malloc2d(double*** array, int n, int m) {
    double* p = reinterpret_cast<double*>(malloc(n * m * sizeof(double)));
    if (!p) return -1;
    (*array) = reinterpret_cast<double**>(malloc(n * sizeof(double*)));
    if (!(*array)) {
        free(p);
        return -1;
    }
    for (int i = 0; i < n; i++)
        (*array)[i] = &(p[i * m]);
    return 0;
}

int free2d(double*** array) {
    free(&((*array)[0][0]));
    free(*array);
    return 0;
}

std::vector<double> Get_max_parallel(double** matrix, int stolb, int stroka) {
    int g;
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    if (ProcNum == 1 || stolb < ProcNum) {
        return Get_max_sequestional( matrix, stolb, stroka);
    }
    std::vector<double> max(ProcRank == 0 ? stolb : stolb / ProcNum);
    int offset = stolb / ProcNum;
    int rem = stolb % ProcNum;
    double* mass = new double[stolb];
    std::map<int, double> map;
    MPI_Status Status;
    double ProcSum = 0.0;
    MPI_Bcast(mass, stolb, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = ProcRank * offset; (i < ProcRank * offset + offset) && (i < stolb); i++) {
        map[i] = findmax(i, matrix, stroka);
    }
    if (ProcRank == 0) {
        for (int i = 0; i < offset; i++) {
            mass[i] = static_cast<double>(map[i]);
        }
        for (int i = 1; i < ProcNum; i++) {
            for (int j = i * offset; (j < i * offset + offset) && (j < stolb); j++) {
                MPI_Recv(&g, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &Status);
                MPI_Recv(&ProcSum, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &Status);
                mass[g] = ProcSum;
            }
        }
    } else {
        for (int i = ProcRank * offset; (i < ProcRank * offset + offset) && (i < stolb); i++) {
            MPI_Send(&i, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
            MPI_Send(&map[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    if (ProcRank == 0) {
        if (rem != 0) {
            for (int j = stolb - rem; j < stolb; j++)
                mass[j] = findmax(j, matrix, stroka);
        }
        for (int i = 0; i < stolb; i++)
            max[i] = mass[i];
    }

    delete[] mass;
    return max;
}

double** Matrix(int stolb, int stroka) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    double** ptrarray = 0;
    malloc2d(&ptrarray, stroka, stolb);
    for (int count_row = 0; count_row < stroka; count_row++)
        for (int count_column = 0; count_column < stolb; count_column++)
            ptrarray[count_row][count_column] = gen() % 10 + 1.0 / static_cast<double>((gen() % 10 + 1.0));
    return ptrarray;
}

std::vector<double> Get_max_sequestional(double** matrix, int stolb, int stroka) {
    std::vector<double> max(stolb);
    for (int i = 0; i < stolb; i++)
       max[i] = findmax(i, matrix, stroka);
    return max;
}
