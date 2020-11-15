// Copyright 2020 Zoreev Mikhail
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <cmath>
#include <limits>

#include "./rectangle_method.h"

TEST(Parallel_Operations_MPI, INTEGRAL_FROM_0_TO_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double parallel_result = integralParallel(static_cast<double(*)(double)>(&cos), 0, 1, 10000);
    if (rank == 0) {
        double sequential_result = integralSeqential(static_cast<double(*)(double)>(&cos), 0, 1, 10000);
        ASSERT_LT(std::fabs(parallel_result - sequential_result), std::numeric_limits<double>::epsilon() * 1000);
    }
}

TEST(Parallel_Operations_MPI, INTEGRAL_FROM_5_TO_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double parallel_result = integralParallel(static_cast<double(*)(double)>(&cos), 5, 0, 10000);
    if (rank == 0) {
        double sequential_result = integralSeqential(static_cast<double(*)(double)>(&cos), 5, 0, 10000);
        ASSERT_LT(std::fabs(parallel_result - sequential_result), std::numeric_limits<double>::epsilon() * 1000);
    }
}

TEST(Parallel_Operations_MPI, INTEGRAL_FROM_0_TO_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double parallel_result = integralParallel(static_cast<double(*)(double)>(&cos), 0, 100, 10000);
    if (rank == 0) {
        double sequential_result = integralSeqential(static_cast<double(*)(double)>(&cos), 0, 100, 10000);
        ASSERT_LT(std::fabs(parallel_result - sequential_result), std::numeric_limits<double>::epsilon() * 1000);
    }
}

TEST(Parallel_Operations_MPI, INTEGRAL_FROM_0_TO_709) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double parallel_result = integralParallel(static_cast<double(*)(double)>(&cos), 0, 709, 10000);
    if (rank == 0) {
        double sequential_result = integralSeqential(static_cast<double(*)(double)>(&cos), 0, 709, 10000);
        ASSERT_LT(std::fabs(parallel_result - sequential_result), std::numeric_limits<double>::epsilon() * 10000);
    }
}

TEST(Parallel_Operations_MPI, INTEGRAL_WITH_LOW_RANGE) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double parallel_result = integralParallel(static_cast<double(*)(double)>(&cos), 1, 1.01, 10000);
    if (rank == 0) {
        double sequential_result = integralSeqential(static_cast<double(*)(double)>(&cos), 1, 1.01, 10000);
        ASSERT_LT(std::fabs(parallel_result - sequential_result), std::numeric_limits<double>::epsilon() * 1000);
    }
}

TEST(Parallel_Operations_MPI, EXCETION_ON_ZERO_COUNT) {
    EXPECT_THROW(integralParallel(static_cast<double(*)(double)>(&cos), 5, 0, 0), std::runtime_error);
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
