// Copyright 2020 Napylov Evgenii
#include <iostream>
#include "./count_diff_char_mpi.h"
#include <mpi.h> //FOR TESTING WITHOUT GTEST !!! DELETE BEFORE PUSH

int main(int argc, char** argv) {
    const char* str1 = "Hello, world!";
    const char* str2 = "hello, human!***";
    std::cout << str1 << std::endl;
    std::cout << str2 << std::endl;
    MPI_Init(NULL, NULL);
    int diff_count = getParallelDifferenceCount(str1, str2);
    MPI_Finalize();
    std::cout << diff_count << std::endl;
    return 0;
}