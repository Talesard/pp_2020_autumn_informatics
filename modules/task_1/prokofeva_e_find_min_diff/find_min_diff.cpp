// Copyright 2020 Prokofeva Elizaveta
#include <mpi.h>
#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include "random"
#include "iostream"
#include "../../../modules/task_1/prokofeva_e_find_min_diff/find_min_diff.h"

std::vector<int> getRandomVector(int size) {
    std::mt19937 random_generator;
    random_generator.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(size);
    for (int i = 0; i < size; i++) { vec[i] = random_generator() % 100; }
    return vec;
}

pair<int, int> getParallelOperations(const vector<int>& vec) {
    int size, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if ((static_cast<int>(vec.size()) <= size) || (size == 1)) {
        if (rank == 0) {
            return getSequentialOperations(vec);
        } else {
            return std::make_pair(0, 0);
        }
    }
    vector<int> local_vector;

    int delta = vec.size() / size + 1;

    if (rank == 0) {
        for (int num_proc = 1; num_proc < size; num_proc++) {
            int start_index = delta * num_proc - 1 * num_proc;

            if (num_proc == size - 1) {
                MPI_Send(&vec[0] + start_index, vec.size() - start_index,
                    MPI_INT, num_proc, 0, MPI_COMM_WORLD);
            } else {
                MPI_Send(&vec[0] + start_index, delta, MPI_INT, num_proc, 0, MPI_COMM_WORLD);
            }
        }
        local_vector = vector<int>(delta);
        local_vector = vector<int>(vec.begin(), vec.begin() + delta);
    } else {
        MPI_Status status;

        if (rank == size - 1) {
            delta = vec.size() - (delta * (size - 1) - 1 * (size - 1));
        }
        local_vector = vector<int>(delta);
        MPI_Recv(&local_vector[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }

    auto result = getSequentialOperations(local_vector);
    int tmp_arr[3] = { result.first, result.second, rank };
    int min_elements[3] = { 0, 0, 0 };

    MPI_Op op;
    MPI_Op_create(reinterpret_cast<MPI_User_function*>(pairCompare), true, &op);
    MPI_Reduce(tmp_arr, min_elements, 3, MPI_INT, op, 0, MPI_COMM_WORLD);
    return std::make_pair(min_elements[0], min_elements[1]);
}

pair<int, int> getSequentialOperations(const vector<int>& vec) {
    if (vec.size() < 2) {
        throw std::runtime_error("Size of vector less than 2!");
    }
    int difference = abs(vec[0] - vec[1]);
    pair<int, int> result = std::make_pair(vec[0], vec[1]);

    for (int i = 1; i < static_cast<int>(vec.size()) - 1; i++) {
        if (abs(vec[i] - vec[i + 1]) < difference) {
            difference = abs(vec[i] - vec[i + 1]);
            result.first = vec[i];
            result.second = vec[i + 1];
        }
    }

    return result;
}

void pairCompare(int* input, int* output, int* lenght, MPI_Datatype* dtype) {
    if (abs(input[0] - input[1]) < abs(output[0] - output[1])) {
        output[0] = input[0];
        output[1] = input[1];
    } else {
        if (abs(input[0] - input[1]) == abs(output[0] - output[1]) && (input[2] < output[2])) {
            output[0] = input[0];
            output[1] = input[1];
        }
    }
}
