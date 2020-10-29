// Copyright 2020 Loganov Andrei
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <ctime>
#include <numeric>
#include "./matrix_sum_stolb.h"

TEST(TEST_PARALEL_MPI, TEST1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int rows = 4;
    int cols = 4;
    int size = rows * cols;
    std::vector<int> matrvec(size, 0);
    std::vector<int> result(cols, 0);
    std::vector<int> result2(cols, 0);
    matrvec = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 };
    result = { 4, 8, 12, 16 };
    result2 = getParallelMaxElements(matrvec, cols, rows);
    if (rank == 0) {
        for (int i = 0; i < cols; i++) {
            ASSERT_EQ(result2[i], result[i]) << "Arrays check_array and res_array differ at index " << i;
        }
    }
}
TEST(TEST_PARALEL_MPI, TEST2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int rows = 5;
    const int cols = 2;
    int size = rows * cols;
    int** matr;
    std::vector<int> matrvec(size, 0);
    std::vector<int> result(cols, 0);
    std::vector<int> result2(cols, 0);
    matr = new int* [rows];
    for (int i = 0; i < rows; i++) {
        matr[i] = new int[cols];
    }
    matrvec = randomMatrix(matr, cols, rows);
    if (rank == 0) {
        result2 = getSequentialMaxStolb(matrvec, cols, rows);
    }
    result = getParallelMaxElements(matrvec, cols, rows);
    if (rank == 0) {
        ASSERT_EQ(result, result2);
    }
}
TEST(TEST_PARALEL_MPI, TEST3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int rows = 2;
    const int cols = 5;
    int size = rows * cols;
    int** matr;
    std::vector<int> matrvec(size, 0);
    std::vector<int> result(cols, 0);
    std::vector<int> result2(cols, 0);
    matr = new int* [rows];
    for (int i = 0; i < rows; i++) {
        matr[i] = new int[cols];
    }
    matrvec = randomMatrix(matr, cols, rows);
    if (rank == 0) {
        result2 = getSequentialMaxStolb(matrvec, cols, rows);
    }
    result = getParallelMaxElements(matrvec, cols, rows);
    if (rank == 0) {
        ASSERT_EQ(result, result2);
    }
}
TEST(TEST_PARALEL_MPI, TEST4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int rows = 134;
    const int cols = 111;
    int size = rows * cols;
    int** matr;
    std::vector<int> matrvec(size, 0);
    std::vector<int> result(cols, 0);
    std::vector<int> result2(cols, 0);
    matr = new int* [rows];
    for (int i = 0; i < rows; i++) {
        matr[i] = new int[cols];
    }
    matrvec = randomMatrix(matr, cols, rows);
    if (rank == 0) {
        result2 = getSequentialMaxStolb(matrvec, cols, rows);
    }
    result = getParallelMaxElements(matrvec, cols, rows);
    if (rank == 0) {
        ASSERT_EQ(result, result2);
    }
}
TEST(TEST_PARALEL_MPI, TEST5) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int rows = 33;
    const int cols = 3;
    int size = rows * cols;
    int** matr;
    std::vector<int> matrvec(size, 0);
    std::vector<int> result(cols, 0);
    std::vector<int> result2(cols, 0);
    matr = new int* [rows];
    for (int i = 0; i < rows; i++) {
        matr[i] = new int[cols];
    }
    matrvec = randomMatrix(matr, cols, rows);
    if (rank == 0) {
        result2 = getSequentialMaxStolb(matrvec, cols, rows);
    }
    result = getParallelMaxElements(matrvec, cols, rows);
    if (rank == 0) {
        ASSERT_EQ(result, result2);
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
