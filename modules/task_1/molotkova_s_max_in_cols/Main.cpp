// Copyright 2020 Molotkova Svetlana
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <vector>
#include "./Max.h"

TEST(Parallel_Operations_MPI, Matrix_6x6) {
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int stolb = 6;
    int stroka = 6;
    double** matrix = Matrix(stolb, stroka);
    std::vector<double> global_vec;
    std::vector<double> result = Get_max_parallel(matrix, stolb, stroka);

    if (ProcRank == 0) {
        std::vector<double> control_result = Get_max_sequestional(matrix, stolb, stroka);
        ASSERT_EQ(control_result, result);
    }
    free2d(&matrix);
}


TEST(Parallel_Operations_MPI, Matrix_1x6) {
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int stolb = 1;
    int stroka = 6;

    double** matrix = Matrix(stolb, stroka);
    std::vector<double> global_vec;
    std::vector<double> result = Get_max_parallel(matrix, stolb, stroka);

    if (ProcRank == 0) {
        std::vector<double> control_result = Get_max_sequestional(matrix, stolb, stroka);
        ASSERT_EQ(control_result, result);
    }
    free2d(&matrix);
}

TEST(Parallel_Operations_MPI, Matrix_3x1) {
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int stolb = 3;
    int stroka = 1;
    double** matrix = Matrix(stolb, stroka);
    std::vector<double> global_vec;
    std::vector<double> result = Get_max_parallel(matrix, stolb, stroka);

    if (ProcRank == 0) {
        std::vector<double> control_result = Get_max_sequestional(matrix, stolb, stroka);
        ASSERT_EQ(control_result, result);
    }

    free2d(&matrix);
}

TEST(Parallel_Operations_MPI, Matrix_19x21) {
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int stolb = 19;
    int stroka = 21;
    double** matrix = Matrix(stolb, stroka);
    std::vector<double> global_vec;
    std::vector<double> result = Get_max_parallel(matrix, stolb, stroka);

    if (ProcRank == 0) {
        std::vector<double> control_result = Get_max_sequestional(matrix, stolb, stroka);
        ASSERT_EQ(control_result, result);
    }

    free2d(&matrix);
}

TEST(Parallel_Operations_MPI, Matrix_14x2) {
    int ProcNum, ProcRank;

    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

    int stolb = 14;
    int stroka = 2;
    double** matrix = Matrix(stolb, stroka);
    std::vector<double> global_vec;
    std::vector<double> result = Get_max_parallel(matrix, stolb, stroka);

    if (ProcRank == 0) {
        std::vector<double> control_result = Get_max_sequestional(matrix, stolb, stroka);
        ASSERT_EQ(control_result, result);
    }

    free2d(&matrix);
}

int main(int argc, char* argv[]) {
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
