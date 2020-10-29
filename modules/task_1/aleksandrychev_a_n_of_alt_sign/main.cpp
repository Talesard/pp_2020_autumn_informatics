// Copyright 2020 Aleksandychev Andrey
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <string>
#include <vector>
#include "./n_of_alt_sign.h"

TEST(Parallel_Count_Sentences_MPI, Rand_Vector_With_Size_50) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 50;
    std::vector<int> global_vector(count);

    if (rank == 0) {
        global_vector = create_rand_vector(count);
    }

    int parallel_res = get_parallel_solution(global_vector);

    if (rank == 0) {
        int single_res = get_num_alter_signs(global_vector);
        ASSERT_EQ(parallel_res, single_res);
    }
}

TEST(Parallel_Count_Sentences_MPI, Rand_Vector_With_Size_500) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 500;
    std::vector<int> global_vector(count);

    if (rank == 0) {
        global_vector = create_rand_vector(count);
    }

    int parallel_res = get_parallel_solution(global_vector);

    if (rank == 0) {
        int single_res = get_num_alter_signs(global_vector);
        ASSERT_EQ(parallel_res, single_res);
    }
}

TEST(Parallel_Count_Sentences_MPI, Rand_Vector_With_Size_50000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 50000;
    std::vector<int> global_vector(count);

    if (rank == 0) {
        global_vector = create_rand_vector(count);
    }

    int parallel_res = get_parallel_solution(global_vector);

    if (rank == 0) {
        int single_res = get_num_alter_signs(global_vector);
        ASSERT_EQ(parallel_res, single_res);
    }
}

TEST(Parallel_Count_Sentences_MPI, Rand_Vector_With_Size_100000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 100000;
    std::vector<int> global_vector(count);

    if (rank == 0) {
        global_vector = create_rand_vector(count);
    }

    int parallel_res = get_parallel_solution(global_vector);

    if (rank == 0) {
        int single_res = get_num_alter_signs(global_vector);
        ASSERT_EQ(parallel_res, single_res);
    }
}

TEST(Parallel_Count_Sentences_MPI, Rand_Vector_With_Size_1000000) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int count = 1000000;
    std::vector<int> global_vector(count);

    if (rank == 0) {
        global_vector = create_rand_vector(count);
    }

    int parallel_res = get_parallel_solution(global_vector);

    if (rank == 0) {
        int single_res = get_num_alter_signs(global_vector);
        ASSERT_EQ(parallel_res, single_res);
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
