// Copyright 2020 Napylov Evgenii
#include <mpi.h>
#include <random>
#include <ctime>
#include "./count_diff_char_mpi.h"


char* getRandomString(int len) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    const char alph[] = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
    char* res = new char[len + 1];
    for (int i = 0; i < len; i++) {
        res[i] = alph[gen() % (strlen(alph) - 1)];
    }
    res[len] = '\0';
    return res;
}

int getSequentialDifferenceCount(const char *str1, const char *str2) {
    int diff_count = 0;
    for (int i = 0; i < strlen(str1); i++) {
        if (str1[i] != str2[i]) {
            diff_count++;
        }
    }
    return diff_count;
}

int getParallelDifferenceCount(const char *str1, const char *str2) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const unsigned len1 = strlen(str1);
    const unsigned len2 = strlen(str2);
    int min_len = 0;
    int len_diff = 0;
    int global_diff = 0;

    if (len1 < len2) {
        min_len = len1;
        len_diff = len2 - len1;
    }
    else {
        min_len = len2;
        len_diff = len1 - len2;
    }

    int rem = min_len % size;
    int delta = min_len / size;

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&str1[0] + rem + proc * delta, delta, MPI_CHAR, proc, 0, MPI_COMM_WORLD);
            MPI_Send(&str2[0] + rem + proc * delta, delta, MPI_CHAR, proc, 1, MPI_COMM_WORLD);
        }
    }

    char* str1_local;
    char* str2_local;
    int local_diff = 0;

    if (rank == 0) {
        str1_local = new char[rem + delta + 1];
        str2_local = new char[rem + delta + 1];
        for (int i = 0; i < rem + delta; i++) { 
            str1_local[i] = str1[i];
            str2_local[i] = str2[i];
        }
        str1_local[rem + delta] = '\0';
        str2_local[rem + delta] = '\0';
    }
    else {
        str1_local = new char[delta + 1];
        str2_local = new char[delta + 1];
        MPI_Status status;
        MPI_Recv(&str1_local[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&str2_local[0], delta, MPI_CHAR, 0, 1, MPI_COMM_WORLD, &status);
        str1_local[delta] = '\0';
        str2_local[delta] = '\0';
    }
    local_diff = getSequentialDifferenceCount(str1_local, str2_local);
    delete[] str1_local;
    delete[] str2_local;
    MPI_Reduce(&local_diff, &global_diff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_diff + len_diff;
}