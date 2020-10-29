// Copyright 2020 Aleksandychev Andrey
#include <mpi.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <random>
#include "../../../../modules/task_1/aleksandrychev_a_n_of_alt_sign/n_of_alt_sign.h"

std::vector<int> create_rand_vector(int size) {
    std::mt19937 random_generator;
    random_generator.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> a;
    for (int i = 0; i < size; i++) { a.push_back(random_generator() % 200 - 100); }
    return a;
}

int get_num_alter_signs(std::vector<int> a) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int res = 0;
    for (int i = 0; i < static_cast<int>(a.size()) - 1; i++) {
        if ((a[i] & (1 << 31)) != (a[i + 1] & (1 << 31))) {
            res++;
        }
    }

    return res;
}

int get_parallel_solution(std::vector<int> a) {
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (static_cast<int>(a.size()) <= size || size == 1) {
        if (rank == 0) {
            return get_num_alter_signs(a);
        } else {
            return -2147483647;
        }
    }

    int delta = a.size() / size;
    int remainder = a.size() % size;

    std::vector<int> locale_vector;
    int loc_res, greate_res = 0;

    if (rank == 0) {
        for (int idx_p = 1; idx_p < size; idx_p++) {
            if (idx_p != size - 1) {
                MPI_Send(&a[0] + delta * idx_p + remainder, delta + 1, MPI_INT, idx_p, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&a[0] + delta * idx_p + remainder, delta, MPI_INT, idx_p, 0, MPI_COMM_WORLD);
            }
        }
    }

    if (rank == 0) {
        locale_vector = std::vector<int>(a.begin(), a.begin() + delta + remainder + 1);
    } else {
        MPI_Status status;

        if (rank != size - 1) {
            locale_vector = std::vector<int>(delta + 1);
        } else {
            locale_vector = std::vector<int>(delta);
        }

        if (rank != size - 1) {
            MPI_Recv(&locale_vector[0], delta + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        } else {
            MPI_Recv(&locale_vector[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        }
    }
    loc_res = get_num_alter_signs(locale_vector);
    MPI_Reduce(&loc_res, &greate_res, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return greate_res;
}
