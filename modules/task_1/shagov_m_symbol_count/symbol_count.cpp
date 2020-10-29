// Copyright 2020 Shagov Maksim
#include <mpi.h>
#include <random>
#include <ctime>
#include <string>
#include "../../../modules/task_1/shagov_m_symbol_count/symbol_count.h"
constexpr char string_for_random[] =
"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890!,.&^$#@*()-_+=#%&'*/:;<>?][~}{|";

std::string createRandomString(int sz) {
    std::string string;
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    for (int i = 0; i < sz; i++) {
        string.push_back(string_for_random[gen() % 92]);
        }
    return string;
}

int calculateSymbolsCountSequental(const std::string& string, char symbol) {
    if (string.size() == 0)
        throw "Invalid size";
    int count = 0;
    for (auto it = string.cbegin(); it != string.cend(); ++it)
        if (*it == symbol)
            count++;
    return count;
}

int getParallelSymbolsCount(const std::string& global_string,
                            int count_size_vector, char symbol) {
    if (global_string.size() == 0)
        throw "Invalid size";
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = count_size_vector / size;
    int remainder = count_size_vector - delta * size;

    if (rank == 0) {
        for (int proc = 1; proc < size; proc++) {
            MPI_Send(&global_string[0] + remainder + proc * delta, delta,
                MPI_CHAR, proc, 0, MPI_COMM_WORLD);
        }
    }
    std::string local_string;
    local_string.resize(delta + (rank == 0 ? remainder : 0));
    if (rank == 0) {
        local_string = std::string(global_string.begin(),
                                   global_string.begin() + delta + remainder);
    } else {
        MPI_Status status;
        MPI_Recv(&local_string[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }

    unsigned int global_sum = 0;
    unsigned int local_sum = calculateSymbolsCountSequental(local_string, symbol);
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum;
}
