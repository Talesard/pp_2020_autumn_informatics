// Copyright 2020 Luckianchencko Ivan
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include<cmath>
#include "./integral.h"

TEST(Paral, Test_func1_1) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func1, 2, 5);
    double answer = get_Integral(func1, 2, 5);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}
TEST(Paral, Test_func1_2) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func1, 1, 10);
    double answer = get_Integral(func1, 1, 10);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}
TEST(Paral, Test_func2_1) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func2, 2, 5);
    double answer = get_Integral(func2, 2, 5);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}
TEST(Paral, Test_func2_2) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func2, 1, 15);
    double answer = get_Integral(func2, 1, 15);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}
TEST(Paral, Test_func3_1) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func3, 2, 5);
    double answer = get_Integral(func3, 2, 5);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}
TEST(Paral, Test_func3_2) {
    int ProcRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double paral_answer = get_Paral_Integral(func3, 2, 15);
    double answer = get_Integral(func3, 2, 15);
    if (ProcRank == 0) {
        ASSERT_DOUBLE_EQ(paral_answer, answer);
    }
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
