// Copyright 2020 Volkova Anastasia
#include <mpi.h>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>
#include "../../../modules/task_1/volkova_a_sents_in_line/sents_mpi.h"

std::vector<char> getRandomLine(int size) {
    std::vector<char> line(size);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    int symbol;
    for (int i = 0; i < size; i++) {
        symbol = gen() % 255;
        line[i] = static_cast<char>(symbol);
    }
    return line;
}

int getSequentialSumSentences(std::vector<char> line) {
    const int  size = line.size();
    int ans = 0;
    for (int i = 0; i < size; i++) {
        if (line[i] == '!' || line[i] == '?' || line[i] == '.') ans++;
    }
    return ans;
}

int getParallelSumSentences(std::vector<char> global_line, int size_line) {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int delta = size_line / size;
    const int tmp = delta * size;
    const int k = size_line - tmp;
    if (rank == 0) {
        for (int process = 1; process < size; process++) {
            MPI_Send(&global_line[0] + k + process * delta, delta, MPI_CHAR, process, 0, MPI_COMM_WORLD);
        }
    }
    int vec_local_size = 0;
    if (rank == 0) {
        vec_local_size = k + delta;
    } else {
        vec_local_size = delta;
    }
    std::vector<char> local_line(vec_local_size);
    if (rank == 0) {
        local_line = std::vector<char>(global_line.begin(), global_line.begin() + k + delta);
    } else {
        MPI_Status status;
        MPI_Recv(&local_line[0], delta, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &status);
    }
    int global_sum_sentences = 0;
    int local_sum_sentences = getSequentialSumSentences(local_line);
    MPI_Reduce(&local_sum_sentences, &global_sum_sentences, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    return global_sum_sentences;
}
