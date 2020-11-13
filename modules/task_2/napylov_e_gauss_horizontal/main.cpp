// Copyright 2020 Napylov Evgenii
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>
#include <vector>
#include "./gauss_horizontal.h"

const double EPSILON = std::numeric_limits<double>::epsilon() * 100;


TEST(Gauss_horizontal_ribbon_MPI, Test_2x3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(2 * 3);
    sys_eq = {
        1, 2,   5,
        3, 4,   7
    };
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, 2, 3);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, 2, 3, res, EPSILON);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        SolveGaussSeq(sys_eq, 2, 3);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
    }
}

TEST(Gauss_horizontal_ribbon_MPI, Test_5x6) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(5 * 6);
    sys_eq = {
         7, 3, 8, 1, 2,     6,
         8, 2, 6, 9, 1,     2,
         2, 4, 5, 1, 7,     1,
         5, 8, 1, 8, 9,     4,
         2, 8, 4, 6, 5,     8
    };
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, 5, 6);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, 5, 6, res, EPSILON);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        SolveGaussSeq(sys_eq, 5, 6);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
