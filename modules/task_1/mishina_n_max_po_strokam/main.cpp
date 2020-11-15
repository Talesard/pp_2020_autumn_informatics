// Copyright 2020 Mishina Nadezhda
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./max_po_strokam.h"

TEST(Parallel_Operations_MPI, Test_3_str_3_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 3;
    int stb = 3;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
TEST(Parallel_Operations_MPI, Test_1000_str_1000_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 1000;
    int stb = 1000;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
TEST(Parallel_Operations_MPI, Test_46_str_115_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 46;
    int stb = 115;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
TEST(Parallel_Operations_MPI, Test_95_str_4_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 95;
    int stb = 4;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
TEST(Parallel_Operations_MPI, Test_10_str_1_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 10;
    int stb = 1;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
TEST(Parallel_Operations_MPI, Test_1_str_10_stb) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int str = 1;
    int stb = 10;
    std::vector<std::vector<int>> check_m(str*stb);
    check_m = random_matrix(str, stb);
    std::vector<int> check_v = convert_to_vector(check_m, str, stb);
    std::vector<int> max(str);
    std::vector<int> check_max(str);
    check_max = max_vector_sequental(check_v, str, stb);
    max = max_vector_parallel(check_v, str, stb);
    if (rank == 0) {
        ASSERT_EQ(check_max, max);
    }
}
int main(int argc, char**argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners = ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
