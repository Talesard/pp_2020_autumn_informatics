// Copyright 2020 Napylov Evgenii
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <limits>
#include <vector>
#include "./gauss_horizontal.h"

const double EPSILON = std::numeric_limits<double>::epsilon() * 100;



TEST(Gauss_horizontal_MPI, Test_2x3) {
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

TEST(Gauss_horizontal_MPI, Test_4x5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(4 * 5);
    sys_eq = {
        5, 1, 7, 8,     5,
        5, 4, 3, 1,     7,
        5, 8, 3, 4,     9,
        1, 6, 4, 9,     3
    };
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, 4, 5);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, 4, 5, res, EPSILON);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        SolveGaussSeq(sys_eq, 4, 5);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
    }
}

TEST(Gauss_horizontal_MPI, Test_5x6) {
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

TEST(Gauss_horizontal_MPI, Test_8x9) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(8 * 9);
    sys_eq = {
        4, 2, 7, 3, 10, 5, 7, 1,    4,
        10, 1, 5, 7, 9, 3, 3, 10,   5,
        9, 5, 10, 7, 7, 3, 3, 1,    9,
        10, 8, 9, 2, 6, 9, 6, 5,    3,
        6, 1, 6, 2, 1, 5, 8, 5,     8,
        1, 8, 6, 2, 9, 8, 6, 4,     3,
        3, 3, 5, 8, 8, 8, 6, 7,     8,
        9, 4, 4, 3, 2, 4, 7, 4,     1
    };
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, 8, 9);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, 8, 9, res, EPSILON);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        SolveGaussSeq(sys_eq, 8, 9);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
    }
}

TEST(Gauss_horizontal_MPI, Test_10x11) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<double> sys_eq(10 * 11);
    sys_eq = {
        4, 3, 10, 8, 8, 1, 7, 10, 10, 1,    10,
        8, 1, 4, 10, 2, 7, 8, 5, 4, 5,      2,
        10, 9, 3, 4, 8, 4, 4, 5, 3, 5,      2,
        2, 4, 5, 5, 3, 6, 6, 6, 2, 3,       2,
        9, 10, 10, 8, 1, 3, 2, 9, 9, 8,     6,
        2, 5, 2, 5, 4, 5, 7, 10, 5, 8,      2,
        4, 6, 6, 7, 5, 3, 3, 8, 8, 5,       8,
        10, 2, 3, 1, 9, 2, 5, 7, 6, 3,      7,
        2, 6, 5, 1, 10, 1, 9, 2, 3, 6,      3,
        3, 4, 6, 4, 4, 4, 2, 9, 1, 6,       5
    };
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, 10, 11);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, 10, 11, res, EPSILON);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        SolveGaussSeq(sys_eq, 10, 11);
        t1 = MPI_Wtime();
        std::cout << "seq_time: " << t1 - t0 << std::endl;
    }
}

TEST(Gauss_horizontal_MPI, Test_performance_180x181) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 180;
    int cols = 181;
    std::vector<double> sys_eq(rows * cols);
    sys_eq = getSpecialMatrix(rows, cols);
    double t0, t1;
    t0 = MPI_Wtime();
    std::vector<double> res = SolveGaussParallel(sys_eq, rows, cols);
    t1 = MPI_Wtime();
    if (rank == 0) {
        bool check = CheckSolution(sys_eq, rows, cols, res, EPSILON * 1e7);
        ASSERT_TRUE(check);
        std::cout << "par_time: " << t1 - t0 << std::endl;
        t0 = MPI_Wtime();
        res = SolveGaussSeq(sys_eq, rows, cols);
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
