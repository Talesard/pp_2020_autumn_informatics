// Copyright 2020 Rustamov Azer
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include "./matrix_cols_sum.h"

TEST(Parallel_Matrix_Cols_Sum, Size_0x0) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 0, cols = 0;
    Matrix matrix;
    if (procRank == 0)
        matrix = RandomMatrix(rows, cols);
    ASSERT_ANY_THROW(ParallelColsSum(matrix, rows, cols));
}

TEST(Parallel_Matrix_Cols_Sum, Size_72x1) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 72, cols = 1;
    Matrix matrix;
    if (procRank == 0)
        matrix = RandomMatrix(rows, cols);
    Matrix sum = ParallelColsSum(matrix, rows, cols);
    if (procRank == 0) {
        Matrix control_sum = SequentialColsSum(matrix, rows, cols);
        for (int i = 0; i < cols; i++) {
            ASSERT_EQ(control_sum[i], sum[i]);
        }
    }
}

TEST(Parallel_Matrix_Cols_Sum, Size_1x100) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 1, cols = 100;
    Matrix matrix;
    if (procRank == 0)
         matrix = RandomMatrix(rows, cols);
    Matrix sum = ParallelColsSum(matrix, rows, cols);
    if (procRank == 0) {
        Matrix control_sum = SequentialColsSum(matrix, rows, cols);
        for (int i = 0; i < cols; i++) {
            ASSERT_EQ(control_sum[i], sum[i]);
        }
    }
}

TEST(Parallel_Matrix_Cols_Sum, Size_22x44) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 22, cols = 44;
    Matrix matrix;
    if (procRank == 0)
        matrix = RandomMatrix(rows, cols);
    Matrix sum = ParallelColsSum(matrix, rows, cols);
    if (procRank == 0) {
        Matrix control_sum = SequentialColsSum(matrix, rows, cols);
        for (int i = 0; i < cols; i++) {
            ASSERT_EQ(control_sum[i], sum[i]);
        }
    }
}

TEST(Parallel_Matrix_Cols_Sum, Size_58x77) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 58, cols = 77;
    Matrix matrix;
    if (procRank == 0)
        matrix = RandomMatrix(rows, cols);
    Matrix sum = ParallelColsSum(matrix, rows, cols);
    if (procRank == 0) {
        Matrix control_sum = SequentialColsSum(matrix, rows, cols);
        for (int i = 0; i < cols; i++) {
            ASSERT_EQ(control_sum[i], sum[i]);
        }
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
