// Copyright 2020 Yurin Stanislav
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_1/yurin_s_sum_of_values_by_matrix_rows/sum_of_values_by_matrix_rows.h"

std::vector<int> getSequentialSumOfMatrixRows(std::vector<int> matrix, int number_of_rows) {
    int matrix_size = matrix.size();
    int number_of_elements_in_one_row;
    std::vector<int> result_vector(number_of_rows);

    if (number_of_rows == 0) {
        result_vector.resize(1, 0);
        return result_vector;
    } else if (matrix_size == 0) {
        result_vector.resize(number_of_rows, 0);
        return result_vector;
    }

    if ((matrix_size % number_of_rows) == 0) {
        number_of_elements_in_one_row  = matrix_size / number_of_rows;
    } else {
        number_of_elements_in_one_row  = (matrix_size / number_of_rows) + 1;
        matrix.resize(number_of_elements_in_one_row * number_of_rows, 0);
    }

    for (int  i = 0; i < number_of_rows; i++) {
        result_vector[i] =
                getSumOfOneRow(std::vector<int>(matrix.begin() + (i * number_of_elements_in_one_row),
                                                (matrix.begin() + (i * number_of_elements_in_one_row)) +
                                                                        number_of_elements_in_one_row));
    }

    return result_vector;
}

std::vector<int> getRandomVector(int sz) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)));
    std::vector<int> vec(sz);
    for (int  i = 0; i < sz; i++) { vec[i] = gen() % 100; }
    return vec;
}

int getSumOfOneRow(std::vector<int> row) {
    const int  sz = row.size();
    int accumulating_elem = 0;

    for (int  i = 0; i < sz; i++) {
        accumulating_elem += row[i];
    }

    return accumulating_elem;
}

std::vector<int> getSumOfOneProc(std::vector<int> row, int number_of_elements_in_one_row) {
    const int  sz = row.size();
    int number_of_rows_in_one_proc = sz / number_of_elements_in_one_row;
    std::vector<int> return_vec_of_one_proc(number_of_rows_in_one_proc);

    for (int  i = 0; i < number_of_rows_in_one_proc; i++) {
        return_vec_of_one_proc[i] =
                getSumOfOneRow(std::vector<int>(row.begin() + (i * number_of_elements_in_one_row),
                                                (row.begin() + (i * number_of_elements_in_one_row)) +
                                                                    number_of_elements_in_one_row));
    }

    return return_vec_of_one_proc;
}

std::vector<int> getParallelSumOfMatrixRows(std::vector<int> matrix,
                                            int number_of_rows) {
    int size, rank;
    int matrix_size = matrix.size();
    int number_of_elements_in_one_row;
    int number_of_rows_for_one_proc_int;
    int number_of_rows_for_one_proc_this;
    int number_of_overload_proc;
    std::vector<int> result_vector(number_of_rows);

    if (number_of_rows == 0) {
        result_vector.resize(1, 0);
        return result_vector;
    } else if (matrix_size == 0) {
        result_vector.resize(number_of_rows, 0);
        return result_vector;
    }

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);



    if ((matrix_size % number_of_rows) == 0) {
        number_of_elements_in_one_row  = matrix_size / number_of_rows;
    } else {
        number_of_elements_in_one_row  = (matrix_size / number_of_rows) + 1;

        matrix.resize(number_of_elements_in_one_row * number_of_rows, 0);
    }

    number_of_rows_for_one_proc_int = number_of_rows / size;
    number_of_overload_proc = number_of_rows % size;
    if (rank < number_of_overload_proc)
        number_of_rows_for_one_proc_this = number_of_rows_for_one_proc_int + 1;
    else
        number_of_rows_for_one_proc_this = number_of_rows_for_one_proc_int;
    std::vector<int> sendcounts(size);
    for (int i = 0; i < size; i++) {
        if (i < number_of_overload_proc)
            sendcounts[i] = (number_of_rows_for_one_proc_int + 1) * number_of_elements_in_one_row;
        else
            sendcounts[i] = number_of_rows_for_one_proc_int * number_of_elements_in_one_row;
    }
    std::vector<int> displs_for_scatterv(size);
    int delta_for_scatterv = 0;
    for (int i = 0; i < size; i++) {
        displs_for_scatterv[i] = delta_for_scatterv;
        delta_for_scatterv += sendcounts[i];
    }

    std::vector<int> rows_for_one_proc(sendcounts[rank]);

    MPI_Scatterv(&matrix[0], &sendcounts[0], &displs_for_scatterv[0], MPI_INT,
                    &rows_for_one_proc[0], sendcounts[rank], MPI_INT,
                    0, MPI_COMM_WORLD);

    std::vector<int> sum_of_one_proc(number_of_rows_for_one_proc_this);
    sum_of_one_proc = getSumOfOneProc(rows_for_one_proc, number_of_elements_in_one_row);

    std::vector<int> recvcounts(size);
    for (int i = 0; i < size; i++) {
        if (i < number_of_overload_proc)
            recvcounts[i] = (number_of_rows_for_one_proc_int + 1);
        else
            recvcounts[i] = number_of_rows_for_one_proc_int;
    }
    std::vector<int> displs_for_gatherv(size);
    int delta_for_gatherv = 0;
    for (int i = 0; i < size; i++) {
        displs_for_gatherv[i] = delta_for_gatherv;
        delta_for_gatherv += recvcounts[i];
    }

    MPI_Gatherv(&sum_of_one_proc[0],
                number_of_rows_for_one_proc_this,
                MPI_INT,
                &result_vector[0],
                &recvcounts[0],
                &displs_for_gatherv[0],
                MPI_INT,
                0, MPI_COMM_WORLD);

    return result_vector;
}
