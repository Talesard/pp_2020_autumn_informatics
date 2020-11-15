// Copyright 2020 Kasyanychev Mikhail
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./kasyanychev_m_min_elem.h"

TEST(Parallel_Operations_MPI, Line_Size_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size = 10;
    std::vector<int> vec(size);
    if (rank == 0) {
        vec = RandomVector(size);
    }
    int parallel = ParallelOperations(vec, size);
    if (rank == 0) {
        int sequential = SequentialOperations(vec);
        ASSERT_EQ(sequential, parallel);
    }
}

TEST(Parallel_Operations_MPI, Line_Size_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size = 100;
    std::vector<int> vec(size);
    if (rank == 0) {
        vec = RandomVector(size);
    }
    int parallel = ParallelOperations(vec, size);
    if (rank == 0) {
        int sequential = SequentialOperations(vec);
        ASSERT_EQ(sequential, parallel);
    }
}

TEST(Parallel_Operations_MPI, Line_Size_1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size = 1000;
    std::vector<int> vec(size);
    if (rank == 0) {
        vec = RandomVector(size);
    }
    int parallel = ParallelOperations(vec, size);
    if (rank == 0) {
        int sequential = SequentialOperations(vec);
        ASSERT_EQ(sequential, parallel);
    }
}

TEST(Parallel_Operations_MPI, Line_Size_10000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size = 10000;
    std::vector<int> vec(size);
    if (rank == 0) {
        vec = RandomVector(size);
    }
    int parallel = ParallelOperations(vec, size);
    if (rank == 0) {
        int sequential = SequentialOperations(vec);
        ASSERT_EQ(sequential, parallel);
    }
}

TEST(Parallel_Operations_MPI, Own_vector) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> vec = { 58, 234, 2014, 23, 845 };
    int parallel = ParallelOperations(vec, vec.size());
    if (rank == 0) {
        int sequential = SequentialOperations(vec);
        ASSERT_EQ(sequential, parallel);
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
