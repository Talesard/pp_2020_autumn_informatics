// Copyright 2020 Kumbrasev Mark
#include <mpi.h>
#include <cstdlib>
#include <string>
#include <iostream>
#include <cmath>
#include "../../../modules/task_1/kumbrasev_m_diff_char_counter/diff_char_counter.h"

std::string get_random_string() {
    static std::mt19937 random;
    std::string str{};
    for (std::size_t k = 0; k < 10000; k++) {
        str += static_cast<char>('a' + random() % 20);
    }
    return str;
}

std::size_t difference_count(seq_policy, const std::string& str_lhs, const std::string& str_rhs) {
    std::size_t counter{};
    for (std::size_t i = 0; i < str_lhs.size(); i++) {
        counter += str_lhs[i] != str_rhs[i];
    }

    return counter;
}

std::size_t difference_count(par_policy, std::string str_lhs, std::string str_rhs) {
    int diff = 0, locdiff = 0, sumdiff = 0, locSize = 0, newSize, ProcNum, rank;
    char *str;

    diff = abs(static_cast<int>(str_lhs.size()) - static_cast<int>(str_rhs.size()));
    str_lhs.size() > str_rhs.size() ? str_lhs.resize(str_rhs.size()) : str_rhs.resize(str_lhs.size());
    int b = str_lhs.size();
    str = new char();
    MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        locSize = (b + ProcNum - 1) / ProcNum;
        newSize = (locSize * ProcNum);
        str_lhs += std::string((newSize - b), '0');
        str_rhs += std::string((newSize - b), '0');
        std::string buf = "";
        str = new char[2 * newSize + 1];

        for (int i = 0; i < newSize; i += locSize)
            buf += str_lhs.substr(i, locSize) + str_rhs.substr(i, locSize);

        for (int i = 0; i < newSize * 2 + 1; ++i) {
            str[i] = buf[i];
        }
    }

    MPI_Bcast(&locSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    char *locstr = reinterpret_cast<char *>(malloc(locSize * 2 * sizeof(char)));
    MPI_Scatter(str, locSize * 2, MPI_CHAR, locstr, locSize * 2, MPI_CHAR, 0, MPI_COMM_WORLD);

    for (int i = 0; i < locSize; i++)
        locdiff += locstr[i] != locstr[i + locSize];
    MPI_Reduce(&locdiff, &sumdiff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    return sumdiff + diff;
}
