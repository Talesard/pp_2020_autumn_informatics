// Copyright 2020 Bulychev Vladislav
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include "./integration_trapez.h"

auto f1 = [](double x) { return sin(x); };
auto f2 = [](double x) { return x - 1 * x + 2 * x * x; };


TEST(Integral_Test, Error_N_less_1) {
    EXPECT_ANY_THROW(getIntegral(0, 1, 0, f2));
}

TEST(Integral_Test, Return_correct_answer) {
    double i = getIntegral(0, 3.14, 3000, f1);
    EXPECT_NEAR(i, 2.0, 0.01);
}

TEST(Parallel_Test, Return_correct_answer) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double i = getParallelIntegral(0, 3.14, 3000, f1);
    if (rank == 0) {
        EXPECT_NEAR(i, 2.0, 0.01);
    }
}

TEST(Parallel_Test, From_0_to_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double parallel_integral = getParallelIntegral(0, 1, 1000, f2);
    if (rank == 0) {
        double integral = getIntegral(0, 1, 1000, f2);
        ASSERT_NEAR(parallel_integral, integral, 0.01);
    }
}

TEST(Parallel_Test, From_0_to_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double parallel_integral = getParallelIntegral(0, 10, 1000, f2);
    if (rank == 0) {
        double integral = getIntegral(0, 10, 1000, f2);
        ASSERT_NEAR(parallel_integral, integral, 0.01);
    }
}

TEST(Parallel_Test, From_0_to_1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    double parallel_integral = getParallelIntegral(0, 1000, 1000, f2);
    if (rank == 0) {
        double integral = getIntegral(0, 1000, 1000, f2);
        ASSERT_NEAR(parallel_integral, integral, 0.01);
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
