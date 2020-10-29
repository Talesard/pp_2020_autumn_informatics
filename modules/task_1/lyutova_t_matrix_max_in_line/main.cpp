// Copyright 2020 Lyutova Tanya
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include"./matrix_max_in_line.h"

TEST(Parallel_Operations_MPI, Test_size4x4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> create_matr;
    const int rows = 4;
    const int columns = 4;

    if (rank == 0) {
        create_matr = getRandomMatrix(rows, columns);
    }

    std::vector<int> global = getMaxInLineParallel(create_matr, rows, columns);

    if (rank == 0) {
        std::vector<int> reference = getMaxInLineSequental(create_matr, rows, columns);
        ASSERT_EQ(reference, global);
    }
}

TEST(Parallel_Operations_MPI, Test_size1x1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> create_matr;
    const int rows = 1;
    const int columns = 1;

    if (rank == 0) {
        create_matr = getRandomMatrix(rows, columns);
    }

    std::vector<int> global = getMaxInLineParallel(create_matr, rows, columns);

    if (rank == 0) {
        std::vector<int> reference = getMaxInLineSequental(create_matr, rows, columns);
        ASSERT_EQ(reference, global);
    }
}

TEST(Parallel_Operations_MPI, Test_size0x0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> create_matr;
    const int rows = 0;
    const int columns = 0;

    if (rank == 0) {
        create_matr = getRandomMatrix(rows, columns);
    }
    ASSERT_ANY_THROW(getMaxInLineParallel(create_matr, rows, columns));
}


TEST(Parallel_Operations_MPI, Test_bigSize) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> create_matr;
    const int rows = 123;
    const int columns = 434;

    if (rank == 0) {
        create_matr = getRandomMatrix(rows, columns);
    }

    std::vector<int> global = getMaxInLineParallel(create_matr, rows, columns);

    if (rank == 0) {
        std::vector<int> reference = getMaxInLineSequental(create_matr, rows, columns);
        ASSERT_EQ(reference, global);
    }
}

TEST(Parallel_Operations_MPI, Test_0xBigsize) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> create_matr;
    const int rows = 0;
    const int columns = 434;

    if (rank == 0) {
        create_matr = getRandomMatrix(rows, columns);
    }

    ASSERT_ANY_THROW(getMaxInLineParallel(create_matr, rows, columns));
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
