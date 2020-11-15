// Copyright 2020 Gogov Vladislav
#include <mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./block_striped_columnwise.h"

TEST(Parallel_Mult_Matrix_By_Vector_MPI, Size_10x20_10) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 10, cols = 20, vecSize = 10;
    Matrix matrix;
    std::vector<double> vec;
    if (procRank == 0) {
        matrix = createRandomMatrix(rows, cols);
        vec = createRandomMatrix(vecSize, 1);
    }
    ASSERT_ANY_THROW(multMatrixByVectorParallel(matrix, rows, cols, vec, vecSize));
}

TEST(Parallel_Mult_Matrix_By_Vector_MPI, Size_10x10_10) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 10, cols = 10, vecSize = 10;
    Matrix matrix;
    std::vector<double> vec;
    if (procRank == 0) {
        matrix = createRandomMatrix(rows, cols);
        vec = createRandomMatrix(vecSize, 1);
    }
    double startTime, endTime;
    if (procRank == 0)
        startTime = MPI_Wtime();
    std::vector<double> resultParallel = multMatrixByVectorParallel(matrix, rows, cols, vec, vecSize);
    if (procRank == 0) {
        endTime = MPI_Wtime();
        std::cout << "Time parallel: " << endTime - startTime << std::endl;
    }
    if (procRank == 0) {
        startTime = MPI_Wtime();
        std::vector<double> resultSequential = multMatrixByVectorSequential(matrix, rows, cols, vec, vecSize);
        endTime = MPI_Wtime();
        std::cout << "Time sequential: " << endTime - startTime << std::endl;
        ASSERT_TRUE(assertVectors(resultParallel, resultSequential));
    }
}

TEST(Parallel_Mult_Matrix_By_Vector_MPI, Size_100x100_100) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 100, cols = 100, vecSize = 100;
    Matrix matrix;
    std::vector<double> vec;
    if (procRank == 0) {
        matrix = createRandomMatrix(rows, cols);
        vec = createRandomMatrix(vecSize, 1);
    }
    double startTime, endTime;
    if (procRank == 0)
        startTime = MPI_Wtime();
    std::vector<double> resultParallel = multMatrixByVectorParallel(matrix, rows, cols, vec, vecSize);
    if (procRank == 0) {
        endTime = MPI_Wtime();
        std::cout << "Time parallel: " << endTime - startTime << std::endl;
    }
    if (procRank == 0) {
        startTime = MPI_Wtime();
        std::vector<double> resultSequential = multMatrixByVectorSequential(matrix, rows, cols, vec, vecSize);
        endTime = MPI_Wtime();
        std::cout << "Time sequential: " << endTime - startTime << std::endl;
        std::cout << "Test status: ";
        ASSERT_TRUE(assertVectors(resultParallel, resultSequential));
    }
}

TEST(Parallel_Mult_Matrix_By_Vector_MPI, Size_460x570_570) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 460, cols = 570, vecSize = 570;
    Matrix matrix;
    std::vector<double> vec;
    if (procRank == 0) {
        matrix = createRandomMatrix(rows, cols);
        vec = createRandomMatrix(vecSize, 1);
    }
    double startTime, endTime;
    if (procRank == 0)
        startTime = MPI_Wtime();
    std::vector<double> resultParallel = multMatrixByVectorParallel(matrix, rows, cols, vec, vecSize);
    if (procRank == 0) {
        endTime = MPI_Wtime();
        std::cout << "Time parallel: " << endTime - startTime << std::endl;
    }
    if (procRank == 0) {
        startTime = MPI_Wtime();
        std::vector<double> resultSequential = multMatrixByVectorSequential(matrix, rows, cols, vec, vecSize);
        endTime = MPI_Wtime();
        std::cout << "Time sequential: " << endTime - startTime << std::endl;
        ASSERT_TRUE(assertVectors(resultParallel, resultSequential));
    }
}

TEST(Parallel_Mult_Matrix_By_Vector_MPI, Size_333x444_444) {
    int procRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int rows = 333, cols = 444, vecSize = 444;
    Matrix matrix;
    std::vector<double> vec;
    if (procRank == 0) {
        matrix = createRandomMatrix(rows, cols);
        vec = createRandomMatrix(vecSize, 1);
    }
    double startTime, endTime;
    if (procRank == 0)
        startTime = MPI_Wtime();
    std::vector<double> resultParallel = multMatrixByVectorParallel(matrix, rows, cols, vec, vecSize);
    if (procRank == 0) {
        endTime = MPI_Wtime();
        std::cout << "Time parallel: " << endTime - startTime << std::endl;
    }
    if (procRank == 0) {
        startTime = MPI_Wtime();
        std::vector<double> resultSequential = multMatrixByVectorSequential(matrix, rows, cols, vec, vecSize);
        endTime = MPI_Wtime();
        std::cout << "Time sequential: " << endTime - startTime << std::endl;
        ASSERT_TRUE(assertVectors(resultParallel, resultSequential));
    }
}

int main(int argc, char* argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners &listeners =
        ::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
