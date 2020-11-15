// Copyright 2020 Yurin Stanislav
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./sum_of_values_by_matrix_rows.h"

TEST(Sum_Of_Values_By_Matrix_Rows, Zero_Number_Of_Rows_And_Zero_Matrix_Size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 0, number_of_rows = 0;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        std::vector<int> result_vector_sequential = getSequentialSumOfMatrixRows(matrix, number_of_rows);
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Sum_Of_Values_By_Matrix_Rows, Zero_Number_Of_Rows_And_Matrix_Size_Is_173) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 173, number_of_rows = 0;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        std::vector<int> result_vector_sequential = getSequentialSumOfMatrixRows(matrix, number_of_rows);
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Sum_Of_Values_By_Matrix_Rows, Number_Of_Rows_Is_61_And_Zero_Matrix_Size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 0, number_of_rows = 61;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        std::vector<int> result_vector_sequential = getSequentialSumOfMatrixRows(matrix, number_of_rows);
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Sum_Of_Values_By_Matrix_Rows, Number_Of_Rows_Is_28_And_Matrix_Size_Is_327) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 327, number_of_rows = 28;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        std::vector<int> result_vector_sequential = getSequentialSumOfMatrixRows(matrix, number_of_rows);
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Sum_Of_Values_By_Matrix_Rows, Number_Of_Rows_More_Than_Matrix_Size) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 9, number_of_rows = 37;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        std::vector<int> result_vector_sequential = getSequentialSumOfMatrixRows(matrix, number_of_rows);
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Sum_Of_Values_By_Matrix_Rows, Result_Vec_Size_Is_Equal_Number_Of_Rows) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int matrix_size = 327, number_of_rows = 28;
    std::vector<int> matrix(matrix_size);

    if (rank == 0) {
        matrix = getRandomVector(matrix_size);
    }

    std::vector<int> result_vector_parallel = getParallelSumOfMatrixRows(matrix, number_of_rows);

    if (rank == 0) {
        ASSERT_EQ((int)result_vector_parallel.size(), number_of_rows);
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
