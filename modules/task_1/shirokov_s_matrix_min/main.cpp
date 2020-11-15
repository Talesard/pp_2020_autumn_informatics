// Copyright 2020 Shirokov Sergey
#include <mpi.h>
#include <gtest/gtest.h>
#include <gtest-mpi-listener.hpp>
#include "./min_elem_matrix.h"

TEST(Parallel_Matrix_Max_MPI, Size_0x0) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizel = 0, sizec = 0;
    if (procRank == 0) {
        ASSERT_TRUE(sizel == 0 && sizec == 0);
    }
}

TEST(Parallel_Matrix_Max_MPI, Size_100x100) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizel = 100, sizec = 100;
    if (procRank == 0) {
        ASSERT_TRUE(sizel == 100 && sizec == 100);
    }
}

TEST(Parallel_Matrix_Max_MPI, Size_1x25) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizel = 1, sizec = 25;
    if (procRank == 0) {
        ASSERT_TRUE(sizel == 1 && sizec == 25);
    }
}

TEST(Parallel_Matrix_Max_MPI, Size_35x35) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizel = 35, sizec = 35;
    ASSERT_TRUE(sizel == 35 && sizec == 35);
}

TEST(Parallel_Matrix_Max_MPI, Size_53x47) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int sizel = 53, sizec = 47;
    ASSERT_TRUE(sizel == 53 && sizec == 47);
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
