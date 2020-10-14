// Copyright 2020 Napylov Evgenii
#include <iostream>
#include "./count_diff_char_mpi.h"
#include <mpi.h> //FOR TESTING WITHOUT GTEST !!! DELETE BEFORE PUSH

int main(int argc, char** argv) {
    const char* str1 = "Hello, world!";
    const char* str2 = "hello, human!***";
    MPI_Init(NULL, NULL);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int diff_count = getParallelDifferenceCount(str1, str2);
    if (rank == 0){
        std::cout << str1 << std::endl;
        std::cout << str2 << std::endl;
        std::cout << diff_count << std::endl;
    }
    MPI_Finalize();
    return 0;
}