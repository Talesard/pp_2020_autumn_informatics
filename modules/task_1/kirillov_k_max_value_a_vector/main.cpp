// Copyright 2020 Kirillov Konstantin
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./max_value_a_vector.h"

TEST(Parallel_Operations_MPI, Test_Max_100) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int size = 100;
    std::vector<int> vec;
    if (procRank == 0) {
        vec = getRandomVector(size);
    }
    int globalMax = getParallelMax(vec, size);
    if (procRank == 0) {
        int seqMax = getSequentialMax(vec);
        ASSERT_EQ(seqMax, globalMax);
    }
}
TEST(Parallel_Operations_MPI, Test_Max_10) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int size = 10;
    std::vector<int> vec;
    if (procRank == 0) {
        vec = getRandomVector(size);
    }
    int globalMax = getParallelMax(vec, size);
    if (procRank == 0) {
        int seqMax = getSequentialMax(vec);
        ASSERT_EQ(seqMax, globalMax);
    }
}
TEST(Parallel_Operations_MPI, Test_Max_300) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int size = 300;
    std::vector<int> vec;
    if (procRank == 0) {
        vec = getRandomVector(size);
    }
    int globalMax = getParallelMax(vec, size);
    if (procRank == 0) {
        int seqMax = getSequentialMax(vec);
        ASSERT_EQ(seqMax, globalMax);
    }
}
TEST(Parallel_Operations_MPI, Test_Max_333) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int size = 500;
    std::vector<int> vec;
    if (procRank == 0) {
        vec = getRandomVector(size);
    }
    int globalMax = getParallelMax(vec, size);
    if (procRank == 0) {
        int seqMax = getSequentialMax(vec);
        ASSERT_EQ(seqMax, globalMax);
    }
}
TEST(Parallel_Operations_MPI, Test_Max_30) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int size = 30;
    std::vector<int> vec;
    if (procRank == 0) {
        vec = getRandomVector(size);
    }
    int globalMax = getParallelMax(vec, size);
    if (procRank == 0) {
        int seqMax = getSequentialMax(vec);
        ASSERT_EQ(seqMax, globalMax);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
