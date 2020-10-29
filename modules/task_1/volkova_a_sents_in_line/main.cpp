// Copyright 2020 Volkova Anastasia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./sents_mpi.h"

TEST(Parallel_MPI, Test_Sentences_Length_0) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 0;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
    }
}

TEST(Parallel_MPI, Test_Sentences_Length_1) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 1;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
    }
}

TEST(Parallel_MPI, Test_Sentences_Length_10) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 10;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
    }
}

TEST(Parallel_MPI, Test_Sentences_Length_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 50;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
    }
}

TEST(Parallel_MPI, Test_Sentences_Length_100) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 100;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
    }
}


TEST(Parallel_MPI, Test_Sentences_Length_1000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<char> global_line;
    const int size_line = 1000;
    if (rank == 0) {
        global_line = getRandomLine(size_line);
    }

    int global_sum_sentences = getParallelSumSentences(global_line, size_line);
    if (rank == 0) {
        int reference_sum_sentences = getSequentialSumSentences(global_line);
        ASSERT_EQ(reference_sum_sentences, global_sum_sentences);
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
