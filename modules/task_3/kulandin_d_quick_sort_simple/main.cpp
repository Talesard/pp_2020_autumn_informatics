// Copyright 2020 Kulandin Denis
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include <algorithm>
#include "./quick_sort.h"

TEST(Parallel_MPI, Test_1) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::vector<int> globalMas;
    const int N = 123;
    if (rank == 0) {
        globalMas = randomVector(N);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    double timeParallel = MPI_Wtime();
    auto parallelRes = parallelQuickSort(globalMas);
    timeParallel = MPI_Wtime() - timeParallel;

    if (rank == 0) {
        double timeSequential = MPI_Wtime();
        std::sort(globalMas.begin(), globalMas.end());
        // quickSort(globalMas, 0, N - 1);
        timeSequential = MPI_Wtime() - timeSequential;
        ASSERT_EQ(parallelRes, globalMas);
        printf("Size of array = %d\n", N);
        printf("Parallel time:              %.10f\n", timeParallel);
        printf("Sequential time:            %.10f\n", timeSequential);
    }
}

TEST(Parallel_MPI, Test_2) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 65431;
    if (rank == 0) {
        globalMas = randomVector(N);
    }

    double timeParallel = MPI_Wtime();
    auto parallelRes = parallelQuickSort(globalMas);
    timeParallel = MPI_Wtime() - timeParallel;

    if (rank == 0) {
        double timeSequential = MPI_Wtime();
        std::sort(globalMas.begin(), globalMas.end());
        timeSequential = MPI_Wtime() - timeSequential;
        ASSERT_EQ(parallelRes, globalMas);
        printf("Size of array = %d\n", N);
        printf("Parallel time:              %.10f\n", timeParallel);
        printf("Sequential time:            %.10f\n", timeSequential);
    }
}

TEST(Parallel_MPI, Test_3) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 1337;
    if (rank == 0) {
        globalMas = randomVector(N);
    }

    double timeParallel = MPI_Wtime();
    auto parallelRes = parallelQuickSort(globalMas);
    timeParallel = MPI_Wtime() - timeParallel;

    if (rank == 0) {
        double timeSequential = MPI_Wtime();
        std::sort(globalMas.begin(), globalMas.end());
        timeSequential = MPI_Wtime() - timeSequential;
        ASSERT_EQ(parallelRes, globalMas);
        printf("Size of array = %d\n", N);
        printf("Parallel time:              %.10f\n", timeParallel);
        printf("Sequential time:            %.10f\n", timeSequential);
    }
}

TEST(Parallel_MPI, Test_4) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::vector<int> globalMas;
    const int N = 1000000;
    if (rank == 0) {
        globalMas = randomVector(N);
    }

    double timeParallel = MPI_Wtime();
    auto parallelRes = parallelQuickSort(globalMas);
    timeParallel = MPI_Wtime() - timeParallel;

    if (rank == 0) {
        double timeSequential = MPI_Wtime();
        std::sort(globalMas.begin(), globalMas.end());
        timeSequential = MPI_Wtime() - timeSequential;
        ASSERT_EQ(parallelRes, globalMas);
        printf("Size of array = %d\n", N);
        printf("Parallel time:              %.10f\n", timeParallel);
        printf("Sequential time:            %.10f\n", timeSequential);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    MPI_Init(&argc, &argv);

    ::testing::AddGlobalTestEnvironment(new GTestMPIListener::MPIEnvironment);
    ::testing::TestEventListeners& listeners =::testing::UnitTest::GetInstance()->listeners();

    listeners.Release(listeners.default_result_printer());
    listeners.Release(listeners.default_xml_generator());

    listeners.Append(new GTestMPIListener::MPIMinimalistPrinter);
    return RUN_ALL_TESTS();
}
