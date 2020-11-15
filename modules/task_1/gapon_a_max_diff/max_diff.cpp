// Copyright 2020 Gapon Andrey

#include <mpi.h>
#include <ctime>
#include <cmath>
#include <vector>
#include <random>
#include <exception>
#include "../../../modules/task_1/gapon_a_max_diff/max_diff.h"

std::vector<int> random_vector(int size) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(nullptr)));
    std::vector<int> vec(size);
    for (int i = 0; i < size; i++)
        vec[i] = gen() % 500;
    return vec;
}

int sequential_search(std::vector<int> vec) {
    const int size = static_cast<int>(vec.size());
    if (size < 2) {
        throw "Error size vector";
    }
    int diff = 0;
    for (int i = 0; i < size - 1; i++) {
        if (abs(vec[i] - vec[i + 1]) > diff)
            diff = abs(vec[i] - vec[i + 1]);
    }
    return diff;
}

int parallel_search(std::vector<int> vec) {
    int vector_size = static_cast<int>(vec.size());
    if (vector_size < 2) {
        throw "Error size vector";
    }
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int delta = (vector_size - 1) / size;
    int remainder = (vector_size - 1) % size;
    std::vector<int> loc_vec;

    if (rank == 0) {
        if (delta) {
            for (int proc = 1; proc < size; proc++) {
                MPI_Send(&vec[0] + proc * delta + remainder, delta + 1, MPI_INT, proc, 0, MPI_COMM_WORLD);
            }
        }
    }

    int diff = 0;
    int diff_max = 0;

    if (rank == 0) {
        if (size != 1 || delta) {
            loc_vec.resize(delta + remainder + 1);
            loc_vec.assign(vec.begin(), vec.begin() + delta + remainder + 1);
            diff_max = sequential_search(loc_vec);
        } else {
            diff_max = sequential_search(vec);
        }
    } else if (delta) {
        loc_vec.resize(delta + 1);
        MPI_Status status;
        MPI_Recv(&loc_vec[0], delta + 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        diff_max = sequential_search(loc_vec);
    }

    MPI_Reduce(&diff_max, &diff, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

    return diff;
}
