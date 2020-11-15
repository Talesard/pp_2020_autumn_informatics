// Copyright 2020 Galkin Vladislav
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./min_values_of_strings.h"

TEST(Min_value, Size_0x0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 0, m = 0;
    std::vector<std::vector<int>> matrix;
    if (rank == 0)
        matrix = createRandomMatrix(n, m);
    ASSERT_ANY_THROW(getMinParallel(matrix));
}

TEST(Min_value, Size_333x287) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 333, m = 287;
    std::vector<std::vector<int>> matrix;
    matrix = createRandomMatrix(n, m);
    std::vector<int> parmin = getMinParallel(matrix);
    if (rank == 0) {
        std::vector<int> seqmin = getMinSequential(matrix);
        ASSERT_EQ(parmin, seqmin);
    }
}

TEST(Min_value, Size_53x1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 53, m = 1;
    std::vector<std::vector<int>> matrix;
    matrix = createRandomMatrix(n, m);
    std::vector<int> parmin = getMinParallel(matrix);
    if (rank == 0) {
        std::vector<int> seqmin = getMinSequential(matrix);
        ASSERT_EQ(parmin, seqmin);
    }
}

TEST(Min_value, Size_1x46) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 1, m = 46;
    std::vector<std::vector<int>> matrix;
    matrix = createRandomMatrix(n, m);
    std::vector<int> parmin = getMinParallel(matrix);
    if (rank == 0) {
        std::vector<int> seqmin = getMinSequential(matrix);
        ASSERT_EQ(parmin, seqmin);
    }
}

TEST(Min_value, Size_104x35) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 104, m = 35;
    std::vector<std::vector<int>> matrix;
    matrix = createRandomMatrix(n, m);
    std::vector<int> parmin = getMinParallel(matrix);
    if (rank == 0) {
        std::vector<int> seqmin = getMinSequential(matrix);
        ASSERT_EQ(parmin, seqmin);
    }
}

TEST(Min_value, Size_62x71) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int n = 62, m = 71;
    std::vector<std::vector<int>> matrix;
    matrix = createRandomMatrix(n, m);
    std::vector<int> parmin = getMinParallel(matrix);
    if (rank == 0) {
        std::vector<int> seqmin = getMinSequential(matrix);
        ASSERT_EQ(parmin, seqmin);
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
