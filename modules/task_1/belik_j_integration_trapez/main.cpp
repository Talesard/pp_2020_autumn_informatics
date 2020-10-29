// Copyright 2020 Belik Julia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include "./IntegrationTrapez.h"

TEST(Parallel_Operations_MPI, Test_Square) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(double)>f = [](double x) {return x * x; };
    double i2 = getParallelIntegration(f, 2, 5, 3);
    if (rank == 0) {
        double i1 = getSequentialIntegration(f, 2, 5, 3);
        ASSERT_NEAR(i1, i2, 1e-2);
    }
}

TEST(Parallel_Operations_MPI, Test_Sin) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(double)>f = [](double x) {return sin(x); };
    double i2 = getParallelIntegration(f, 0, 1, 10);
    if (rank == 0) {
        double i1 = getSequentialIntegration(f, 0, 1, 10);
        EXPECT_NEAR(i1, i2, 1e-2);
    }
}

TEST(Parallel_Operations_MPI, Test_Polyn) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(double)>f = [](double x) {return x * x - 2 * x - 5; };
    double i2 = getParallelIntegration(f, -2, 4, 6);
    if (rank == 0) {
        double i1 = getSequentialIntegration(f, -2, 4, 6);
        EXPECT_NEAR(i1, i2, 1e-2);
    }
}

TEST(Parallel_Operations_MPI, Test_Cube) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(double)>f = [](double x) { return x * x * x; };
    double i2 = getParallelIntegration(f, -100, 100, 100);
    if (rank == 0) {
        double i1 = getSequentialIntegration(f, -100, 100, 100);
        EXPECT_NEAR(i1, i2, 1e-2);
    }
}

TEST(Parallel_Operations_MPI, Test_Cube_One_Interval) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::function<double(double)>f = [](double x) { return x * x * x; };
    double i2 = getParallelIntegration(f, -100, 100, 1);
    if (rank == 0) {
        double i1 = getSequentialIntegration(f, -100, 100, 1);
        EXPECT_NEAR(i1, i2, 1e-2);
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
