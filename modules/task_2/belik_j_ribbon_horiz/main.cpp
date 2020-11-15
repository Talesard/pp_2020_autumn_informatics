// Copyright 2020 Belik Julia
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <mpi.h>
#include <algorithm>
#include <functional>
#include <vector>
#include "./RibbonHoriz.h"

TEST(Parallel_Operations_MPI, 6x6) {
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int rows = 6;
    int cols = 6;
    std::vector<double> matrix;
    if (ProcRank == 0) {
        matrix = Matrix(rows, cols);
    }
    std::vector<double> tmp_vec = Vector(rows);
    std::vector<double> v1 = Get_parallel(matrix, rows, cols, tmp_vec);
    if (ProcRank == 0) {
        std::vector<double> v2 = Get_sequestional(matrix, rows, cols, tmp_vec);
        ASSERT_EQ(v1, v2);
    }
}

TEST(Parallel_Operations_MPI, 1x7) {
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int rows = 1;
    int cols = 7;
    std::vector<double> matrix;
    if (ProcRank == 0) {
        matrix = Matrix(rows, cols);
    }
    std::vector<double> tmp_vec = Vector(rows);
    std::vector<double> v1 = Get_parallel(matrix, rows, cols, tmp_vec);
    if (ProcRank == 0) {
        std::vector<double> v2 = Get_sequestional(matrix, rows, cols, tmp_vec);
        ASSERT_EQ(v1, v2);
    }
}

TEST(Parallel_Operations_MPI, 8x1) {
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int rows = 8;
    int cols = 1;
    std::vector<double> matrix;
    if (ProcRank == 0) {
        matrix = Matrix(rows, cols);
    }
    std::vector<double> tmp_vec = Vector(rows);
    std::vector<double> v1 = Get_parallel(matrix, rows, cols, tmp_vec);
    if (ProcRank == 0) {
        std::vector<double> v2 = Get_sequestional(matrix, rows, cols, tmp_vec);
        ASSERT_EQ(v1, v2);
    }
}

TEST(Parallel_Operations_MPI, 3x3) {
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int rows = 3;
    int cols = 3;
    std::vector<double> matrix;
    if (ProcRank == 0) {
        matrix = Matrix(rows, cols);
    }
    std::vector<double> tmp_vec = Vector(rows);
    std::vector<double> v1 = Get_parallel(matrix, rows, cols, tmp_vec);
    if (ProcRank == 0) {
        std::vector<double> v2 = Get_sequestional(matrix, rows, cols, tmp_vec);
        ASSERT_EQ(v1, v2);
    }
}

TEST(Parallel_Operations_MPI, 16x16) {
    int ProcNum, ProcRank;
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    int rows = 16;
    int cols = 16;
    std::vector<double> matrix;
    if (ProcRank == 0) {
        matrix = Matrix(rows, cols);
    }
    std::vector<double> tmp_vec = Vector(rows);
    std::vector<double> v1 = Get_parallel(matrix, rows, cols, tmp_vec);
    if (ProcRank == 0) {
        std::vector<double> v2 = Get_sequestional(matrix, rows, cols, tmp_vec);
        ASSERT_EQ(v1, v2);
    }
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
