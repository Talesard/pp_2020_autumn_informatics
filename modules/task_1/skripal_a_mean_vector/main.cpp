// Copyright 2020 Skripal Andrey
#include<mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./mean_vector.h"


TEST(Parallel_Operations, Test_Mean_10) {
    int rank;
    double globalMean = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 10;
    std::vector<int>vec;
    if (rank == 0) {
        vec = createRandomVector(size);
    }
    double globalMean1 = getMeanVector(vec, size);
    if (rank == 0) {
        globalMean = getSequentialMean(vec);
        ASSERT_DOUBLE_EQ(globalMean, globalMean1);
    }
}

TEST(Parallel_Operations, Test_Mean_23) {
    int rank;
    double globalMean = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 23;
    std::vector<int>vec;
    if (rank == 0) {
        vec = createRandomVector(size);
    }
    double globalMean1 = getMeanVector(vec, size);
    if (rank == 0) {
        globalMean = getSequentialMean(vec);
        ASSERT_DOUBLE_EQ(globalMean, globalMean1);
    }
}

TEST(Parallel_Operations, Test_Mean_55) {
    int rank;
    double globalMean = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 55;
    std::vector<int>vec;
    if (rank == 0) {
        vec = createRandomVector(size);
    }
    double globalMean1 = getMeanVector(vec, size);
    if (rank == 0) {
        globalMean = getSequentialMean(vec);
        ASSERT_DOUBLE_EQ(globalMean, globalMean1);
    }
}

TEST(Parallel_Operations, Test_Mean_80) {
    int rank;
    double globalMean = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 80;
    std::vector<int>vec;
    if (rank == 0) {
        vec = createRandomVector(size);
    }
    double globalMean1 = getMeanVector(vec, size);
    if (rank == 0) {
        globalMean = getSequentialMean(vec);
        ASSERT_DOUBLE_EQ(globalMean, globalMean1);
    }
}

TEST(Parallel_Operations, Test_Mean_101) {
    int rank;
    double globalMean = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 101;
    std::vector<int>vec;
    if (rank == 0) {
        vec = createRandomVector(size);
    }
    double globalMean1 = getMeanVector(vec, size);
    if (rank == 0) {
        globalMean = getSequentialMean(vec);
        ASSERT_DOUBLE_EQ(globalMean, globalMean1);
    }
}



int main(int argc, char* argv[]) {
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
