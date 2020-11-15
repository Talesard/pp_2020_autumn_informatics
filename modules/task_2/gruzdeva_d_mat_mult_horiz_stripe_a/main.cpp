// Copyright 2020 Gruzdeva Diana
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include "./mat_mult_horiz_stripe_a.h"

TEST(Parallel_Operations_MPI, MATRICES_3x3) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 3, aCols = 3, bRows = 3, bCols = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
    }
}

TEST(Parallel_Operations_MPI, MATRICES_10x2_AND_2x7) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 10, aCols = 2, bRows = 2, bCols = 7;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
    }
}

TEST(Parallel_Operations_MPI, MATRICES_50x25_AND_25x47) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 50, aCols = 25, bRows = 25, bCols = 47;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
    }
}


TEST(Parallel_Operations_MPI, MATRICES_100x100) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 100, aCols = 100, bRows = 100, bCols = 100;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
    }
}

TEST(Parallel_Operations_MPI, MATRICES_100x50_AND_50x12) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 100, aCols = 50, bRows = 50, bCols = 12;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
    }
}

TEST(Parallel_Operations_MPI, MATRICES_240x240) {
    int rank;
    std::vector<double> a;
    std::vector<double> b;
    int aRows = 240, aCols = 240, bRows = 240, bCols = 240;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    a = getRandomMatrix(aRows, aCols, time(0));
    b = getRandomMatrix(bRows, bCols, time(0) + 1);
    double start, end;
    if (rank == 0) {
        start = MPI_Wtime();
    }
    std::vector<double> parallel_result = getParalleMultiplication(a, b, aRows, aCols, bRows, bCols);
    if (rank == 0) {
        end = MPI_Wtime();
        std::cout << "Time for parallel is " << end - start << std::endl;
        start = MPI_Wtime();
        std::vector<double> sequential_result = getSequentialMultiplication(a, b, aRows, aCols, bCols);
        end = MPI_Wtime();
        std::cout << "Time for sequential is " << end - start << std::endl;
        ASSERT_EQ(sequential_result, parallel_result);
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
