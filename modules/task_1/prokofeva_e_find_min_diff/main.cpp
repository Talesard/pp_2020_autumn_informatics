// Copyright 2020 Prokofeva Elizaveta
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "./find_min_diff.h"

TEST(Parallel_Operations_MPI, test_getRandomVector_wrong_size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    const int count_size_vector = -100;

    if (rank == 0) {
        ASSERT_ANY_THROW(getRandomVector(count_size_vector));
    }
}

TEST(Parallel_Operations_MPI, test_getSequentialOperations_works) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<int> global_vector;

    const int count_size_vector = 100;

    if (rank == 0) {
        global_vector = getRandomVector(count_size_vector);
    }

    if (rank == 0) {
        ASSERT_NO_THROW(getSequentialOperations(global_vector));
    }
}

TEST(Parallel_Operations_MPI, test_getSequentialOperations_works_correct) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        vector<int> global_vector = { 2, 5, 4 };
        auto result = getSequentialOperations(global_vector);

        EXPECT_EQ(result.first, 5);
        EXPECT_EQ(result.second, 4);
    }
}

TEST(Parallel_Operations_MPI, test_getParallelOperations_works) {
    std::vector<int> global_vec(5);

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = getRandomVector(5);
    }
    ASSERT_NO_THROW(getParallelOperations(global_vec));
}

TEST(Parallel_Operations_MPI, test_getParallelOperations_and_getSequentialOperations_compare) {
    std::vector<int> global_vec(10);

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = getRandomVector(10);
    }

    auto parallel_result = getParallelOperations(global_vec);

    if (ProcRank == 0) {
        auto sequential_result = getSequentialOperations(global_vec);
        ASSERT_EQ(parallel_result, sequential_result);
    }
}

TEST(Parallel_Operations_MPI, test_getParallelOperations_and_getSequentialOperations_compare_four_size) {
    std::vector<int> global_vec(4);

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = getRandomVector(4);
    }

    auto parallel_result = getParallelOperations(global_vec);

    if (ProcRank == 0) {
        auto sequential_result = getSequentialOperations(global_vec);
        ASSERT_EQ(parallel_result, sequential_result);
    }
}
TEST(Parallel_Operations_MPI, test_getParallelOperations_and_getSequentialOperations_compare_average_size) {
    std::vector<int> global_vec(196);

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = getRandomVector(196);
    }

    auto parallel_result = getParallelOperations(global_vec);

    if (ProcRank == 0) {
        auto sequential_result = getSequentialOperations(global_vec);
        ASSERT_EQ(parallel_result, sequential_result);
    }
}
TEST(Parallel_Operations_MPI, test_getParallelOperations_and_getSequentialOperations_compare_big_size) {
    std::vector<int> global_vec(10655);

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = getRandomVector(10655);
    }

    auto parallel_result = getParallelOperations(global_vec);

    if (ProcRank == 0) {
        auto sequential_result = getSequentialOperations(global_vec);
        ASSERT_EQ(parallel_result, sequential_result);
    }
}
TEST(Parallel_Operations_MPI, test_getParallelOperations_and_getSequentialOperations_compare_works_correct) {
    std::vector<int> global_vec = { 10, 2, 3, 4, 6, 5 };

    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    if (ProcRank == 0) {
        global_vec = { 10, 2, 3, 4, 6, 5 };
    }

    auto parallel_result = getParallelOperations(global_vec);

    if (ProcRank == 0) {
        auto sequential_result = getSequentialOperations(global_vec);
        ASSERT_EQ(parallel_result, sequential_result);
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
