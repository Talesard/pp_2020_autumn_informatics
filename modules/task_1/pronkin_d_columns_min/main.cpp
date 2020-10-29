// Copyright 2020 Pronkin Dmitry
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./columns_min.h"

TEST(Parallel_Operations_MPI, Matrix_5x5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 5;
    int cols = 5;
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

TEST(Parallel_Operations_MPI, Matrix_100x100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 100;
    int cols = 100;
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

TEST(Parallel_Operations_MPI, Matrix_57x28) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 57;
    int cols = 28;
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

TEST(Parallel_Operations_MPI, Matrix_51x1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 51;
    int cols = 1;
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

TEST(Parallel_Operations_MPI, Matrix_1x51) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 1;
    int cols = 51;
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
