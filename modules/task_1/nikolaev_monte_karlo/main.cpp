// Copyright 2020 Nikolaev Denis
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include<mpi.h>
#include<cmath>
#include "./Integral.h"

TEST(PARAL_INTEG, INTEG_f1_1) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    double* y;
    x = new double[SIZE];
    y = new double[SIZE];
    getOperations(f1, x, y, &top1, 0, 2, &Spar1);
    double paral_res = calculating_the_paral_integral(f1, x, y, &Spar1);
    if (ProcRank == 0) {
        double res = calculating_the_integral(f1, x, y, &Spar1);
        ASSERT_GE(paral_res, res);
    }
}

TEST(PARAL_INTEG, INTEG_f1_2) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    x = new double[SIZE];
    double* y;
    y = new double[SIZE];
    getOperations(f1, x, y, &top1, 1, 3, &Spar1);
    double paral_res = calculating_the_paral_integral(f1, x, y, &Spar1);
    if (ProcRank == 0) {
        double res = calculating_the_integral(f1, x, y, &Spar1);
        ASSERT_GE(paral_res, res);
    }
}

TEST(PARAL_INTEG, INTEG_f2_1) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    x = new double[SIZE];
    double* y;
    y = new double[SIZE];
    getOperations(f2, x, y, &top1, 0, 2, &Spar1);
    double paral_res = calculating_the_paral_integral(f2, x, y, &Spar1);
    if (ProcRank == 0) {
        double res = calculating_the_integral(f2, x, y, &Spar1);
        ASSERT_GE(paral_res, res);
    }
}

TEST(PARAL_INTEG, INTEG_f2_2) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    x = new double[SIZE];
    double* y;
    y = new double[SIZE];
    getOperations(f2, x, y, &top1, 2, 4, &Spar1);
    double paral_res = calculating_the_paral_integral(f2, x, y, &Spar1);
    if (ProcRank == 0) {
        double res = calculating_the_integral(f2, x, y, &Spar1);
        ASSERT_GE(paral_res, res);
    }
}

TEST(PARAL_INTEG, INTEG_f3_1) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    x = new double[SIZE];
    double* y;
    y = new double[SIZE];
    getOperations(f3, x, y, &top1, 0, 2, &Spar1);
    double paral_res = calculating_the_paral_integral(f3, x, y, &Spar1);
    double res = calculating_the_integral(f3, x, y, &Spar1);
    if (ProcRank == 0) {
        ASSERT_GE(paral_res, res);
    }
}

TEST(PARAL_INTEG, INTEG_f3_2) {
    int ProcRank;
    int SIZE = 20000;
    MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
    double top1 = 0;
    double Spar1 = 0;
    double* x;
    x = new double[SIZE];
    double* y;
    y = new double[SIZE];
    getOperations(f3, x, y, &top1, 1, 3, &Spar1);
    double paral_res = calculating_the_paral_integral(f3, x, y, &Spar1);
    double res = calculating_the_integral(f3, x, y, &Spar1);
    if (ProcRank == 0) {
        ASSERT_GE(paral_res, res);
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
