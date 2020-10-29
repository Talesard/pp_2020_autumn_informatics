// Copyright 2020 Kulandin Denis
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./matrix_sum.h"

TEST(Parallel_MPI, Test_Sum_10_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 10;
    const int M = 1;
    const int globalSize = N * M;
    if (rank == 0) {
        globalMas = randomVector(globalSize);
    }

    int64_t globalSum = getParallelSum(globalMas, globalSize);
    if (rank == 0) {
        int64_t seqSum = getSequentialSum(globalMas);
        ASSERT_EQ(seqSum, globalSum);
    }
}

TEST(Parallel_MPI, Test_Sum_128_512) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 128;
    const int M = 512;
    const int globalSize = N * M;
    if (rank == 0) {
        globalMas = randomVector(globalSize);
    }

    int64_t globalSum = getParallelSum(globalMas, globalSize);
    if (rank == 0) {
        int64_t seqSum = getSequentialSum(globalMas);
        ASSERT_EQ(seqSum, globalSum);
    }
}

TEST(Parallel_MPI, Test_Sum_1_5000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 1;
    const int M = 5000;
    const int globalSize = N * M;
    if (rank == 0) {
        globalMas = randomVector(globalSize);
    }

    int64_t globalSum = getParallelSum(globalMas, globalSize);
    if (rank == 0) {
        int64_t seqSum = getSequentialSum(globalMas);
        ASSERT_EQ(seqSum, globalSum);
    }
}

TEST(Parallel_MPI, Test_Sum_500_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 500;
    const int M = 1;
    const int globalSize = N * M;
    if (rank == 0) {
        globalMas = randomVector(globalSize);
    }

    int64_t globalSum = getParallelSum(globalMas, globalSize);
    if (rank == 0) {
        int64_t seqSum = getSequentialSum(globalMas);
        ASSERT_EQ(seqSum, globalSum);
    }
}

TEST(Parallel_MPI, Test_Sum_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 0;
    const int M = 1337;
    const int globalSize = N * M;
    if (rank == 0) {
        globalMas = randomVector(globalSize);
    }

    int64_t globalSum = getParallelSum(globalMas, globalSize);
    if (rank == 0) {
        int64_t seqSum = getSequentialSum(globalMas);
        ASSERT_EQ(seqSum, globalSum);
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
