// Copyright 2020 Gapon Andrey

#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./max_diff.h"

TEST(MAX_DIFF_IN_VECTOR_MPI, INCORRECT_SIZE) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vector_size = 1;
    std::vector<int> vec(vector_size);

    if (rank == 0)
        vec = random_vector(vector_size);

    ASSERT_ANY_THROW(parallel_search(vec));
}

TEST(MAX_DIFF_IN_VECTOR_MPI, SIZE_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vector_size = 3;
    std::vector<int> vec(vector_size);

    if (rank == 0)
        vec = random_vector(vector_size);

    int parallel_diff = parallel_search(vec);

    if (rank == 0) {
        int sequential_diff = sequential_search(vec);
        EXPECT_EQ(sequential_diff, parallel_diff);
    }
}

TEST(MAX_DIFF_IN_VECTOR_MPI, SIZE_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vector_size = 50;
    std::vector<int> vec(vector_size);

    if (rank == 0)
        vec = random_vector(vector_size);

    int parallel_diff = parallel_search(vec);

    if (rank == 0) {
        int sequential_diff = sequential_search(vec);
        EXPECT_EQ(sequential_diff, parallel_diff);
    }
}

TEST(MAX_DIFF_IN_VECTOR_MPI, SIZE_500) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vector_size = 500;
    std::vector<int> vec(vector_size);

    if (rank == 0)
        vec = random_vector(vector_size);

    int parallel_diff = parallel_search(vec);

    if (rank == 0) {
        int sequential_diff = sequential_search(vec);
        EXPECT_EQ(sequential_diff, parallel_diff);
    }
}

TEST(MAX_DIFF_IN_VECTOR_MPI, SIZE_5000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int vector_size = 5000;
    std::vector<int> vec(vector_size);

    if (rank == 0)
        vec = random_vector(vector_size);

    int parallel_diff = parallel_search(vec);

    if (rank == 0) {
        int sequential_diff = sequential_search(vec);
        ASSERT_EQ(sequential_diff, parallel_diff);
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
