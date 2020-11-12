// Copyright 2020 Napylov Evgenii
#include <iostream>
#include <limits>
#include <vector>
#include <mpi.h>
#include "./gauss_horizontal.h"

const double EPSILON = std::numeric_limits<double>::epsilon();

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(5*6);
    sys_eq = {
         7, 3, 8, 1, 2,     6,
         8, 2, 6, 9, 1,     2,
         2, 0, 5, 1, 7,     1,
         5, 8, 1, 8, 9,     4,
         2, 8, 4, 6, 5,     8
    };
    std::vector<double> res(5);
    res = SolveGaussParallel(sys_eq, 5, 6);
    std::cout << "res: "; print_vec(res);
    MPI_Finalize();
    return 0;
}
