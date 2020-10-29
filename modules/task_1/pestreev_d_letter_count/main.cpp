// Copyright 2020 Pestreev Daniil
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <string>
#include "./letter_count.h"

TEST(Parallel_Operations_MPI, empty_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string global_str;
    const int size_str = 0;
    if (rank == 0) {
        global_str = getRandomString(size_str);
    }

    int64_t global_sum = LetterInLineParallel(global_str, size_str);
    if (rank == 0) {
        int64_t expexted = 0;
        ASSERT_EQ(expexted, global_sum);
    }
}

TEST(Parallel_Operations_MPI, 14_character_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string global_str;
    const int size_str = 14;
    if (rank == 0) {
        global_str = "3453g453f44d53";
    }

    int64_t global_sum = LetterInLineParallel(global_str, size_str);
    if (rank == 0) {
        int64_t expexted = 3;
        ASSERT_EQ(expexted, global_sum);
    }
}

TEST(Parallel_Operations_MPI, 100_random_character_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string global_str;
    const int size_str = 100;
    if (rank == 0) {
        global_str = getRandomString(size_str);
    }

    int64_t global_sum = LetterInLineParallel(global_str, size_str);
    if (rank == 0) {
        int64_t expexted = LetterInLineSequential(global_str);
        ASSERT_EQ(expexted, global_sum);
    }
}

TEST(Parallel_Operations_MPI, 1000_random_character_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string global_str;
    const int size_str = 1000;
    if (rank == 0) {
        global_str = getRandomString(size_str);
    }

    int64_t global_sum = LetterInLineParallel(global_str, size_str);
    if (rank == 0) {
        int64_t expexted = LetterInLineSequential(global_str);
        ASSERT_EQ(expexted, global_sum);
    }
}

TEST(Parallel_Operations_MPI, 1000_same_character_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int size_str = 1000;
    std::string global_str(size_str, '\n');
    const char G = 'G';
    for (int i = 0; i < size_str; i++) {
        global_str[i] = G;
    }
    int64_t global_sum = LetterInLineParallel(global_str, size_str);

    if (rank == 0) {
        int64_t expexted = LetterInLineSequential(global_str);
        ASSERT_EQ(expexted, global_sum);
    }
}

TEST(Parallel_Operations_MPI, no_letter_character_string) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::string global_str;
    const int size_str = 20;
    if (rank == 0) {
        global_str = "1234567890!@#$%^&*()";
    }

    int64_t global_sum = LetterInLineParallel(global_str, size_str);
    if (rank == 0) {
        int64_t expexted = 0;
        ASSERT_EQ(expexted, global_sum);
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
