// Copyright 2020 Alibekov Murad
#include <mpi.h>
#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/alibekov_m_count_sentences/count_sentences.h"

std::string getRandomString(int sz) {
    if (sz <= 0)
        throw "Invalid size";
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::string str = "";
    for (int i = 0; i < sz; i++) { str += alphabet[static_cast<int>(gen()) % 57]; }
    return str;
}

int getSentencesCountSequential(std::string str) {
    int count = 0;
    if (str.size() == 0)
        return count;
    for (auto elem : str)
        if ((elem == '.') || (elem == '!') || (elem == '?'))
            count++;
    return count;
}

int getSentencesCountParallel(std::string global_str, int count_chars) {
    int proc_count, proc_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &proc_count);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    const int delta = count_chars / proc_count;
    const int remain = count_chars - proc_count * delta;
    if (proc_rank == 0) {
        for (int proc = 1; proc < proc_count; proc++) {
            MPI_Send(global_str.data() + remain + proc * delta, delta, MPI_CHAR, proc, 0, MPI_COMM_WORLD);
        }
    }
    std::string local_str = "";
    if (proc_rank == 0) {
        local_str = std::string(global_str.cbegin(),
                                global_str.cbegin() + delta + remain);
    } else {
        std::vector<char> part_line(delta);
        MPI_Status status;
        MPI_Recv(&part_line[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
        local_str = std::string(part_line.cbegin(), part_line.cend());
    }

    int global_count = 0;
    int local_count = getSentencesCountSequential(local_str);
    MPI_Reduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_count;
}
