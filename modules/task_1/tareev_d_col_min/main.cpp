// Copyright 2020 Tareev Daniil
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./col_min.h"

TEST(Parallel_Operations_MPI, Matrix_2x2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 2;
    int cols = 2;
    std::vector<int> global_vec;
    if (rank == 0) {
        global_vec = transposeMatrix(randomMatrix(cols, rows), cols, rows);
    }
    std::vector<int> result = getParallelOperations(global_vec, rows, cols);
    if (rank == 0) {
        std::vector<int> control_result = getSequentialOperations(global_vec, rows, cols);
        ASSERT_EQ(control_result, result);
    }
}

TEST(Parallel_Operations_MPI, Matrix_10x10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 10;
    int cols = 10;
    std::vector<int> global_vec;
    if (rank == 0) {
        global_vec = transposeMatrix(randomMatrix(cols, rows), cols, rows);
    }
    std::vector<int> result = getParallelOperations(global_vec, rows, cols);
    if (rank == 0) {
        std::vector<int> control_result = getSequentialOperations(global_vec, rows, cols);
        ASSERT_EQ(control_result, result);
    }
}

TEST(Parallel_Operations_MPI, Matrix_75x31) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 75;
    int cols = 31;
    std::vector<int> global_vec;
    if (rank == 0) {
        global_vec = transposeMatrix(randomMatrix(cols, rows), cols, rows);
    }
    std::vector<int> result = getParallelOperations(global_vec, rows, cols);
    if (rank == 0) {
        std::vector<int> control_result = getSequentialOperations(global_vec, rows, cols);
        ASSERT_EQ(control_result, result);
    }
}

TEST(Parallel_Operations_MPI, Matrix_67x44) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 67;
    int cols = 44;
    std::vector<int> global_vec;
    if (rank == 0) {
        global_vec = transposeMatrix(randomMatrix(cols, rows), cols, rows);
    }
    std::vector<int> result = getParallelOperations(global_vec, rows, cols);
    if (rank == 0) {
        std::vector<int> control_result = getSequentialOperations(global_vec, rows, cols);
        ASSERT_EQ(control_result, result);
    }
}

TEST(Parallel_Operations_MPI, Matrix_4x90) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 4;
    int cols = 90;
    std::vector<int> global_vec;
    if (rank == 0) {
        global_vec = transposeMatrix(randomMatrix(cols, rows), cols, rows);
    }
    std::vector<int> result = getParallelOperations(global_vec, rows, cols);
    if (rank == 0) {
        std::vector<int> control_result = getSequentialOperations(global_vec, rows, cols);
        ASSERT_EQ(control_result, result);
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
