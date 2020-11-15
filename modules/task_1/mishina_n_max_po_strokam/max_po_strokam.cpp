// Copyright 2020 Mishina Nadezhda
#include <mpi.h>
#include <vector>
#include <random>
#include<algorithm>
#include <ctime>
#include <limits>
#include "../../../modules/task_1/mishina_n_max_po_strokam/max_po_strokam.h"


std::vector<std::vector<int>> random_matrix(int str, int stb) {
    std::vector<std::vector<int>> matr(str*stb);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < str; i++) {
        for (int j = 0; j < stb; j++) {
            matr[i].push_back(gen());
        }
    }
    return matr;
}
std::vector<int> convert_to_vector(std::vector<std::vector<int>> v, int str, int stb) {
    std::vector<int> res(stb*str);
    for (int i = 0; i < str; i++) {
        for (int j = 0; j < stb; j++) {
            res[i * stb + j] = v[i][j];
        }
    }
    return res;
}
std::vector<int> max_vector_sequental(std::vector<int> v, int str, int stb) {
    std::vector<int> res(str);
    for (int i = 0; i < str; i++) {
        res[i] = *std::max_element(v.begin() + stb*i, v.begin() + stb*(i + 1));
    }
    return res;
}
std::vector<int> max_vector_parallel(std::vector<int> v, int str, int stb) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rem = str % size;
    int delta = str / size;
    int local_count = 0;
    int local_str = 0;
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            MPI_Send(v.data() + i * delta * stb + rem * stb, delta*stb, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        local_count = delta * stb + rem * stb;
        local_str = delta + rem;
    } else {
        local_count = delta * stb;
        local_str = delta;
    }
    std::vector<int> local_v(local_count);
    std::vector<int> local_max(local_str);
    if (rank == 0) {
        local_v = std::vector<int>(v.begin(), v.begin() + local_count);
    } else {
        MPI_Status status;
        MPI_Recv(local_v.data(), local_count, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    }
    if (local_v.size() > 0) {
        local_max = max_vector_sequental(local_v, local_str, stb);
    }
    std::vector<int> max(str);
    if (rank == 0) {
        for (int i = 0; i < local_str; i++) {
            max[i] = local_max[i];
        }
        for (int i = 1; i < size; i++) {
            MPI_Status status;
            MPI_Recv(max.data() + rem + delta * i, delta, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
        }
    } else {
        MPI_Send(local_max.data(), delta, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    return max;
}
