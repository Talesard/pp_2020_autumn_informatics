// Copyright 2020 Galindo Buitrago Julio Cesar
#include<mpi.h>
#include <gtest-mpi-listener.hpp>
#include <gtest/gtest.h>
#include <vector>
#include "./average_vector.h"


TEST(Parallel_Operations, Test_Mean_12) {
    int rango;
    double PromedioGeneral = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int size = 12;
    std::vector<int>vectorial;
    if (rango == 0) {
        vectorial = VectorAleatorio(size);
    }
    double PromedioGeneraal = VectorMedia(vectorial, size);
    if (rango == 0) {
        PromedioGeneral = IntermedioSigma(vectorial);
        ASSERT_DOUBLE_EQ(PromedioGeneral, PromedioGeneraal);
    }
}

TEST(Parallel_Operations, Test_Mean_37) {
    int rango;
    double PromedioGeneral = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int size = 37;
    std::vector<int>vectorial;
    if (rango == 0) {
        vectorial = VectorAleatorio(size);
    }
    double PromedioGeneraal = VectorMedia(vectorial, size);
    if (rango == 0) {
        PromedioGeneral = IntermedioSigma(vectorial);
        ASSERT_DOUBLE_EQ(PromedioGeneral, PromedioGeneraal);
    }
}

TEST(Parallel_Operations, Test_Mean_40) {
    int rango;
    double PromedioGeneral = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int size = 40;
    std::vector<int>vectorial;
    if (rango == 0) {
        vectorial = VectorAleatorio(size);
    }
    double PromedioGeneraal = VectorMedia(vectorial, size);
    if (rango == 0) {
        PromedioGeneral = IntermedioSigma(vectorial);
        ASSERT_DOUBLE_EQ(PromedioGeneral, PromedioGeneraal);
    }
}

TEST(Parallel_Operations, Test_Mean_69) {
    int rango;
    double PromedioGeneral = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int size = 69;
    std::vector<int>vectorial;
    if (rango == 0) {
        vectorial = VectorAleatorio(size);
    }
    double PromedioGeneraal = VectorMedia(vectorial, size);
    if (rango == 0) {
        PromedioGeneral = IntermedioSigma(vectorial);
        ASSERT_DOUBLE_EQ(PromedioGeneral, PromedioGeneraal);
    }
}

TEST(Parallel_Operations, Test_Mean_91) {
    int rango;
    double PromedioGeneral = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rango);
    int size = 91;
    std::vector<int>vectorial;
    if (rango == 0) {
        vectorial = VectorAleatorio(size);
    }
    double PromedioGeneraal = VectorMedia(vectorial, size);
    if (rango == 0) {
        PromedioGeneral = IntermedioSigma(vectorial);
        ASSERT_DOUBLE_EQ(PromedioGeneral, PromedioGeneraal);
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
