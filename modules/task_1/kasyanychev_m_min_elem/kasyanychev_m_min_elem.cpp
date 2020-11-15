// Copyright 2020 Kasyanychev Mikhail
#include <mpi.h>
#include <cstdlib>
#include <ctime>
#include <vector>
#include "../../../modules/task_1/kasyanychev_m_min_elem/kasyanychev_m_min_elem.h"

std::vector<int> RandomVector(int n) {
    if (n < 0)
        throw "Wrong size";
    srand(time(0));
    std::vector<int> res(n);
    for (int i = 0; i < n; i++) {
        res.push_back(std::rand());
    }
    return res;
}

int SequentialOperations(std::vector<int> x) {
    int min = x[0];
    for (size_t i = 1; i < x.size(); i++) {
        if (x[i] < min)
            min = x[i];
    }
    return min;
}

int ParallelOperations(std::vector<int> vec, int len) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = len / size;
    const int epsilon = len - delta * size;
    if (rank == 0) {
        for (int i = 1; i < size; i++)
            MPI_Send(vec.data() + epsilon + delta * i, delta, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
    std::vector<int> lcl;
    if (rank == 0) {
        std::vector<int> tmp(vec.cbegin(), vec.cbegin() + epsilon + delta);
        lcl = tmp;
    } else {
        std::vector<int> tmp(delta);
        MPI_Status status;
        MPI_Recv(&tmp[0], delta, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        std::vector<int> tmp1(tmp.cbegin(), tmp.cend());
        lcl = tmp1;
    }
    int fir = SequentialOperations(lcl);
    int comm = 0;
    MPI_Reduce(&fir, &comm, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    return comm;
}
