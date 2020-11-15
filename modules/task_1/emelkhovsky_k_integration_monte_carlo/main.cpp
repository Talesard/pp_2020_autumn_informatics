// Copyright 2020 Emelkhovsky Ekaterina
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./integration_monte_carlo.h"


TEST(Parallel_Operations_MPI, Test1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a, b, N;
    a = 1;
    b = 10;
    N = 10000;
    double dop = static_cast<double>(b - a) / N;
    double parallel_result = ParallelFunc(N, a, b) * dop;
    if (rank == 0) {
        double easy_result = EasyFunc(N, a, b);
        ASSERT_NEAR(parallel_result, easy_result, 1);
    }
}

TEST(Parallel_Operations_MPI, Test2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a, b, N;
    a = 1;
    b = 100;
    N = 1000;
    double dop = static_cast<double>(b - a) / N;
    double parallel_result = ParallelFunc(N, a, b) * dop;
    if (rank == 0) {
        double easy_result = EasyFunc(N, a, b);
        ASSERT_NEAR(parallel_result, easy_result, 500);
    }
}

TEST(Parallel_Operations_MPI, Test3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a, b, N;
    a = -3;
    b = -1;
    N = 1000;
    double dop = static_cast<double>(b - a) / N;
    double parallel_result = ParallelFunc(N, a, b) * dop;
    if (rank == 0) {
        double easy_result = EasyFunc(N, a, b);
        ASSERT_NEAR(parallel_result, easy_result, 2);
    }
}

TEST(Parallel_Operations_MPI, Test4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a, b, N;
    a = 1;
    b = 10;
    N = 100;
    double dop = static_cast<double>(b - a) / N;
    double parallel_result = ParallelFunc(N, a, b) * dop;
    if (rank == 0) {
        double easy_result = EasyFunc(N, a, b);
        ASSERT_NEAR(parallel_result, easy_result, 10);
    }
}

TEST(Parallel_Operations_MPI, Test5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int a, b, N;
    a = -1;
    b = 10;
    N = 1000;
    double dop = static_cast<double>(b - a) / N;
    double parallel_result = ParallelFunc(N, a, b) * dop;
    if (rank == 0) {
        double easy_result = EasyFunc(N, a, b);
        ASSERT_NEAR(parallel_result, easy_result, 10);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);
    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =
            ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());
    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
