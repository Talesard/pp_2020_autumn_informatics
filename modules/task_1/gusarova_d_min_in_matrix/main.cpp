// Copyright 2020 Gysarova Daria
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./min_in_matr.h"

TEST(Parallel_Operations_MPI, Test_3x3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 3;
    int st = 3;
    std::vector<std::vector<int>> B(str*st);
    B = RandMatr(str, st);
    std::vector<int> A = MatrInArr(B, str, st);
    int min;
    int real_min = MinInArr(A);
    min = ParMin(A);
    if (rank == 0) {
        ASSERT_EQ(real_min, min);
    }
}
TEST(Parallel_Operations_MPI, Test_30x20) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 30;
    int st = 20;
    std::vector<std::vector<int>> B(str*st);
    B = RandMatr(str, st);
    std::vector<int> A = MatrInArr(B, str, st);
    int min;
    int real_min = MinInArr(A);
    min = ParMin(A);
    if (rank == 0) {
        ASSERT_EQ(real_min, min);
    }
}
TEST(Parallel_Operations_MPI, Test_20x30) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 20;
    int st = 30;
    std::vector<std::vector<int>> B(str*st);
    B = RandMatr(str, st);
    std::vector<int> A = MatrInArr(B, str, st);
    int min;
    int real_min = MinInArr(A);
    min = ParMin(A);
    if (rank == 0) {
        ASSERT_EQ(real_min, min);
    }
}
TEST(Parallel_Operations_MPI, Test_13x1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 13;
    int st = 1;
    std::vector<std::vector<int>> B(str*st);
    B = RandMatr(str, st);
    std::vector<int> A = MatrInArr(B, str, st);
    int min;
    int real_min = MinInArr(A);
    min = ParMin(A);
    if (rank == 0) {
        ASSERT_EQ(real_min, min);
    }
}
TEST(Parallel_Operations_MPI, Test_100x100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 100;
    int st = 100;
    std::vector<std::vector<int>> B(str*st);
    B = RandMatr(str, st);
    std::vector<int> A = MatrInArr(B, str, st);
    int min;
    int real_min = MinInArr(A);
    min = ParMin(A);
    if (rank == 0) {
        ASSERT_EQ(real_min, min);
    }
}
int main(int argc, char**argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
