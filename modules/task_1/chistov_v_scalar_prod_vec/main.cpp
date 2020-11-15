// Copyright 2020 Chistov Vladimir

#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./scalar_prod_vec.h"

TEST(Parallel_Count_Sentences_MPI, Vec_Size_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 0;
    std::vector<int> glob_a(size);
    std::vector<int> glob_b(size);

    if (rank == 0) {
        glob_a = GenRandVec(size);
        glob_b = GenRandVec(size);
    }

    int par_sum = CalcSumParallel(glob_a, glob_b);

    if (rank == 0) {
        int seq_sum = CalcSumSequential(glob_a, glob_b);
        ASSERT_EQ(par_sum, seq_sum);
    }
}

TEST(Parallel_Count_Sentences_MPI, Vec_Size_15) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size = 15;
    std::vector<int> glob_a(size);
    std::vector<int> glob_b(size);
    if (rank == 0) {
        glob_a = GenRandVec(size);
        glob_b = GenRandVec(size);
    }

    int par_sum = CalcSumParallel(glob_a, glob_b);

    if (rank == 0) {
        int seq_sum = CalcSumSequential(glob_a, glob_b);
        ASSERT_EQ(par_sum, seq_sum);
    }
}

TEST(Parallel_Count_Sentences_MPI, Vec_Size_25) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 25;
    std::vector<int> glob_a(size);
    std::vector<int> glob_b(size);

    if (rank == 0) {
        glob_a = GenRandVec(size);
        glob_b = GenRandVec(size);
    }

    int par_sum = CalcSumParallel(glob_a, glob_b);

    if (rank == 0) {
        int seq_sum = CalcSumSequential(glob_a, glob_b);
        ASSERT_EQ(par_sum, seq_sum);
    }
}

TEST(Parallel_Count_Sentences_MPI, Vec_Size_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 50;
    std::vector<int> glob_a(size);
    std::vector<int> glob_b(size);

    if (rank == 0) {
        glob_a = GenRandVec(size);
        glob_b = GenRandVec(size);
    }

    int par_sum = CalcSumParallel(glob_a, glob_b);

    if (rank == 0) {
        int seq_sum = CalcSumSequential(glob_a, glob_b);
        ASSERT_EQ(par_sum, seq_sum);
    }
}

TEST(Parallel_Count_Sentences_MPI, Vec_Size_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = 100;
    std::vector<int> glob_a(size);
    std::vector<int> glob_b(size);

    if (rank == 0) {
        glob_a = GenRandVec(size);
        glob_b = GenRandVec(size);
    }

    int par_sum = CalcSumParallel(glob_a, glob_b);

    if (rank == 0) {
        int seq_sum = CalcSumSequential(glob_a, glob_b);
        ASSERT_EQ(par_sum, seq_sum);
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
