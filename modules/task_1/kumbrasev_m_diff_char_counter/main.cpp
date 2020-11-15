// Copyright 2020 Kumbrasev Mark
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>
#include "./diff_char_counter.h"

TEST(Parallel_Operations_MPI, empty_strings) {
    int rank;
    std::string str1 = "";
    std::string str2 = "";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int seq_res = difference_count(seq_policy{}, str1, str2);
    int par_res = difference_count(par_policy{}, str1, str2);
    if (rank == 0) {
        ASSERT_EQ(par_res, seq_res);
    }
}

TEST(Parallel_Operations_MPI, one_empty_string) {
    int rank;
    int exp = 5;
    std::string str1 = "";
    std::string str2 = "abcde";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int par_res = difference_count(par_policy{}, str1, str2);
    if (rank == 0) {
        ASSERT_EQ(par_res, exp);
    }
}

TEST(Parallel_Operations_MPI, equal_strings) {
    int rank;
    std::string str1 = "abcde";
    std::string str2 = "abcde";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int seq_res = difference_count(seq_policy{}, str1, str2);
    int par_res = difference_count(par_policy{}, str1, str2);
    if (rank == 0) {
        ASSERT_EQ(par_res, seq_res);
    }
}

TEST(Parallel_Operations_MPI, different_strings) {
    int rank;
    std::string str1 = "abcde";
    std::string str2 = "qwert";
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int seq_res = difference_count(seq_policy{}, str1, str2);
    int par_res = difference_count(par_policy{}, str1, str2);
    if (rank == 0) {
        ASSERT_EQ(par_res, seq_res);
    }
}

TEST(Parallel_Operations_MPI, random_strings) {
    int rank;
    std::string str1 = get_random_string();
    std::string str2 = get_random_string();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int seq_res = difference_count(seq_policy{}, str1, str2);
    int par_res = difference_count(par_policy{}, str1, str2);
    if (rank == 0) {
        ASSERT_EQ(par_res, seq_res);
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
