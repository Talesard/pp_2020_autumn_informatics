// Copyright 2020 Alibekov Murad
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <string>
#include "./count_sentences.h"

TEST(Parallel_Count_Sentences_MPI, Line_With_Size_50) {
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    const int count_chars = 50;
    std::string global_str = "";

    if (proc_rank == 0) {
        global_str = getRandomString(count_chars);
    }

    int global_count = getSentencesCountParallel(global_str, count_chars);

    if (proc_rank == 0) {
        int reference_count = getSentencesCountSequential(global_str);
        ASSERT_EQ(reference_count, global_count);
    }
}

TEST(Parallel_Count_Sentences_MPI, Line_With_Size_500) {
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    const int count_chars = 500;
    std::string global_str = "";

    if (proc_rank == 0) {
        global_str = getRandomString(count_chars);
    }

    int global_count = getSentencesCountParallel(global_str, count_chars);

    if (proc_rank == 0) {
        int reference_count = getSentencesCountSequential(global_str);
        ASSERT_EQ(reference_count, global_count);
    }
}

TEST(Parallel_Count_Sentences_MPI, Line_With_Size_5000) {
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    const int count_chars = 5000;
    std::string global_str = "";

    if (proc_rank == 0) {
        global_str = getRandomString(count_chars);
    }

    int global_count = getSentencesCountParallel(global_str, count_chars);

    if (proc_rank == 0) {
        int reference_count = getSentencesCountSequential(global_str);
        ASSERT_EQ(reference_count, global_count);
    }
}

TEST(Parallel_Count_Sentences_MPI, Line_With_Size_50000) {
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    const int count_chars = 50000;
    std::string global_str = "";

    if (proc_rank == 0) {
        global_str = getRandomString(count_chars);
    }

    int global_count = getSentencesCountParallel(global_str, count_chars);

    if (proc_rank == 0) {
        int reference_count = getSentencesCountSequential(global_str);
        ASSERT_EQ(reference_count, global_count);
    }
}

TEST(Parallel_Count_Sentences_MPI, Line_With_My_Text) {
    int proc_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);
    std::string global_str = "Hello, world!";

    int global_count = getSentencesCountParallel(global_str, global_str.size());

    if (proc_rank == 0) {
        int reference_count = getSentencesCountSequential(global_str);
        ASSERT_EQ(reference_count, global_count);
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
