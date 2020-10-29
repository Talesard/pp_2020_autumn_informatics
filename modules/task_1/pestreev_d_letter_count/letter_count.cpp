// Copyright 2020 Pestreev Daniil
#include <mpi.h>
#include <time.h>
#include <random>
#include <string>
#include "../../../modules/task_1/pestreev_d_letter_count/letter_count.h"

std::string getRandomString(int size) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_int_distribution<int> distribution(0, 127);
    std::string str(size, '\0');
    for (int i = 0; i < size; i++) {
        str[i] = static_cast<char>(distribution(gen));
    }
    return str;
}

int64_t LetterInLineSequential(std::string str) {
    int count = 0;
    for (size_t i = 0; i < str.size(); i++) {
        if ((str[i] >= 'A' && str[i] <= 'Z') || (str[i] >= 'a' && str[i] <= 'z')) count++;
    }
    return count;
}

int64_t LetterInLineParallel(std::string global_str, int size_str) {
    if ( size_str == 0 ) { return 0; }
    int64_t globalcount = 0;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = size_str / size;
    int remain = size_str - delta*size;
    if (rank == 0) {
        for (int process = 1; process < size; process++) {
            MPI_Send(&global_str[0] + process * delta + remain, delta, MPI_CHAR, process, 0, MPI_COMM_WORLD);
        }
    }

    std::string local_str(delta, '\0');
    if (rank == 0) {
        local_str.resize(delta + remain);
        for (int i = 0; i < delta + remain; i++) {
            local_str[i] = global_str[i];
        }
    } else {
        MPI_Status status;
        MPI_Recv(&local_str[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
    int64_t local_sum = LetterInLineSequential(local_str);
    MPI_Reduce(&local_sum, &globalcount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return globalcount;
}
