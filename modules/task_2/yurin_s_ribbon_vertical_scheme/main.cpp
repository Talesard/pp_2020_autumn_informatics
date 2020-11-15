// Copyright 2020 Yurin Stanislav
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <iostream>
#include "./ribbon_vertical_scheme.h"

TEST(Ribbon_Vertical_Scheme, Zero_Number_Of_Rows) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_columns;
    int number_of_matrix_rows = 0;
    if (rank == 0) {
        std::mt19937 gen;
        gen.seed(get_random_time());
        number_of_matrix_columns = gen() % 100;
    }
    MPI_Bcast(&number_of_matrix_columns, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Ribbon_Vertical_Scheme, Zero_Number_Of_Columns) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_rows;
    int number_of_matrix_columns = 0;
    if (rank == 0) {
        std::mt19937 gen;
        gen.seed(get_random_time());
        number_of_matrix_rows = gen() % 100;
    }
    MPI_Bcast(&number_of_matrix_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Ribbon_Vertical_Scheme, Zero_Number_Of_Rows_And_Zero_Number_Of_Columns) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_rows = 0;
    int number_of_matrix_columns = 0;

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Ribbon_Vertical_Scheme, Number_Of_Rows_Is_31_And_Number_Of_Columns_Is_7) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_rows, number_of_matrix_columns;
    number_of_matrix_rows = 43;
    number_of_matrix_columns = 17;

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Ribbon_Vertical_Scheme, Number_Of_Rows_Is_1_And_Number_Of_Columns_Is_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_rows, number_of_matrix_columns;
    number_of_matrix_rows = 1;
    number_of_matrix_columns = 1;

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
    }
}

TEST(Ribbon_Vertical_Scheme, Random_Number_Of_Rows_And_Random_Number_Of_Columns) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int number_of_matrix_rows, number_of_matrix_columns;
    if (rank == 0) {
        std::mt19937 gen1, gen2;
        gen1.seed(get_random_time());
        number_of_matrix_rows = gen1() % 100;

        gen2.seed(get_random_time());
        number_of_matrix_columns = gen2() % 100;
    }
    MPI_Bcast(&number_of_matrix_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&number_of_matrix_columns, 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> matrix(number_of_matrix_rows, number_of_matrix_columns);
    std::vector<int> factor_vector(number_of_matrix_columns);

    if (rank == 0) {
        matrix = getRandomMatrix(number_of_matrix_rows, number_of_matrix_columns);
        factor_vector = getRandomVector(number_of_matrix_columns);
    }

    std::vector<int> result_vector_parallel =
        getParallelResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);

    if (rank == 0) {
        std::vector<int> result_vector_sequential =
            getSequentialResultVector(matrix, number_of_matrix_rows, number_of_matrix_columns, factor_vector);;
        ASSERT_EQ(result_vector_parallel, result_vector_sequential);
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
