// Copyright 2020 Yurin Stanislav
#include <mpi.h>
#include <vector>
#include <string>
#include <random>
#include <ctime>
#include <algorithm>
#include "../../../modules/task_2/yurin_s_ribbon_vertical_scheme/ribbon_vertical_scheme.h"

unsigned int time_increase = static_cast<unsigned int>(time(0));

unsigned int get_random_time() {
    return time_increase++;
}

std::vector<int> getRandomVector(int sz) {
    std::mt19937 gen;
    gen.seed(get_random_time());

    std::vector<int> vec(sz);
    for (int  i = 0; i < sz; i++)
        vec[i] = gen() % 100;

    return vec;
}

std::vector<int> getRandomMatrix(int number_of_rows, int number_of_columns) {
    std::mt19937 gen;
    gen.seed(static_cast<unsigned int>(time(0)) + 1);

    int sz = number_of_rows * number_of_columns;
    std::vector<int> vec(sz);
    for (int  i = 0; i < sz; i++)
        vec[i] = gen() % 100;

    return vec;
}

std::vector<int> getSequentialResultVector(std::vector<int> matrix, int number_of_matrix_rows,
                                            int number_of_matrix_columns, std::vector<int> factor_vector) {
    // !!!Если нулевая матрица, то возвращаем нулевой вектор
    if (number_of_matrix_rows == 0 || matrix.size() == 0)
        return std::vector<int>(0);
    // Если нулевая матрица, то возвращаем нулевой вектор!!!

    for (int i = 0; i < number_of_matrix_rows; ++i)
        for (int j = 0; j < number_of_matrix_columns; ++j)
            matrix[j + (i * number_of_matrix_columns)] *= factor_vector[j];

    std::vector<int> result_vector(number_of_matrix_rows);

    for (int i = 0; i < number_of_matrix_rows; ++i)
        for (int j = 0; j < number_of_matrix_columns; ++j)
            result_vector[i] += matrix[j + (i * number_of_matrix_columns)];

    return result_vector;
}

std::vector<int> getParallelResultVector(std::vector<int> matrix, int number_of_matrix_rows,
                                            int number_of_matrix_columns, std::vector<int> factor_vector) {
    // !!!Если нулевая матрица, то возвращаем нулевой вектор
    if (number_of_matrix_rows == 0 || number_of_matrix_columns == 0 || matrix.size() == 0)
        return std::vector<int>(0);
    // Если нулевая матрица, то возвращаем нулевой вектор!!!

    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // !!!Массивы смещений и кол-ва процессов для каждого процесса (для overload и не overload процессов)
    int number_of_overload_proc = number_of_matrix_columns % size;
    int number_of_columns_for_one_proc_int = number_of_matrix_columns / size;
    std::vector<int> displs(size);
    std::vector<int> number_of_columns_for_one_proc(size);

    int delta = 0;
    for (int i = 0; i < size; i++) {
        displs[i] = delta;
        if (i < number_of_overload_proc) {
            delta += number_of_columns_for_one_proc_int + 1;

            number_of_columns_for_one_proc[i] = number_of_columns_for_one_proc_int + 1;
        } else {
            delta += number_of_columns_for_one_proc_int;

            number_of_columns_for_one_proc[i] = number_of_columns_for_one_proc_int;
        }
    }
    // Массивы смещений и кол-ва процессов для каждого процесса (для overload и не overload процессов)!!!

    // !!!Создание новых типов для передачи столбцов (для overload и не overload процессов разные типы)
    MPI_Datatype VerticalRibbonType;

    MPI_Type_vector(number_of_matrix_rows, number_of_columns_for_one_proc_int,
                    number_of_matrix_columns, MPI_INT, &VerticalRibbonType);
    MPI_Type_commit(&VerticalRibbonType);

    MPI_Datatype VerticalRibbonForOverloadProcsType;

    MPI_Type_vector(number_of_matrix_rows, number_of_columns_for_one_proc_int + 1,
                    number_of_matrix_columns, MPI_INT, &VerticalRibbonForOverloadProcsType);
    MPI_Type_commit(&VerticalRibbonForOverloadProcsType);
    // Создание новых типов для передачи столбцов (для overload и не overload процессов разные типы)!!!

    // !!!Передаем каждому процессу свою ленту
    int matrix_size = number_of_matrix_rows * number_of_matrix_columns;
    std::vector<int> ribbon_for_one_proc(matrix_size);

    if (rank == 0) {
        // Заполняем ленту для всех процессов кроме 0
        for (int i = 1; i < size; i++)
            if (number_of_columns_for_one_proc[i] != 0) {
                if (i < number_of_overload_proc)
                    MPI_Send(&matrix[displs[i]], 1, VerticalRibbonForOverloadProcsType, i, 0, MPI_COMM_WORLD);
                else
                    MPI_Send(&matrix[displs[i]], 1, VerticalRibbonType, i, 0, MPI_COMM_WORLD);
            }

        // Заполняем ленту для 0 процесса
        for (int i = 0; i < number_of_matrix_rows; i++)
            for (int j = 0; j < number_of_columns_for_one_proc[0]; ++j)
                ribbon_for_one_proc[j + (i * number_of_matrix_columns)] = matrix[j + (i * number_of_matrix_columns)];
    }

    if (rank != 0) {
        MPI_Status status;
        if (number_of_columns_for_one_proc[rank] != 0) {
            if (rank < number_of_overload_proc)
                MPI_Recv(&ribbon_for_one_proc[displs[rank]], 1, VerticalRibbonForOverloadProcsType,
                            0, 0, MPI_COMM_WORLD, &status);
            else
                MPI_Recv(&ribbon_for_one_proc[displs[rank]], 1, VerticalRibbonType, 0, 0, MPI_COMM_WORLD, &status);
        }
    }
    // Передаем каждому процессу свою ленту!!!

    // !!!Передаем элементы factor_vector'а каждому процессу
    std::vector<int> factor_elements(number_of_columns_for_one_proc[rank] + 1);  // + 1 для процессов, раб. в холостую

    MPI_Scatterv(&factor_vector[0], &number_of_columns_for_one_proc[0], &displs[0], MPI_INT,
                    &factor_elements[0], number_of_columns_for_one_proc[rank], MPI_INT, 0, MPI_COMM_WORLD);
    // Передаем элементы factor_vector'а каждому процессам!!!

    // !!!Умножение каждого столбца ленты на свой элемент factor_vector'а
    for (int i = 0; i < number_of_matrix_rows; ++i)
        for (int j = 0; j < number_of_columns_for_one_proc[rank]; ++j)
            ribbon_for_one_proc[displs[rank] + (j + (i * number_of_matrix_columns))] *=
                factor_elements[j];
    // Умножение каждого столбца ленты на свой элемент factor_vector'а!!!

    // !!!Схлопывание лент внутри каждого процесса
    std::vector<int> result_vector_of_one_proc(number_of_matrix_rows);
    for (int i = 0; i < number_of_matrix_rows; ++i)
        for (int j = 0; j < number_of_columns_for_one_proc[rank]; ++j)
            result_vector_of_one_proc[i] +=
                ribbon_for_one_proc[displs[rank] + (j + (i * number_of_matrix_columns))];
    // Схлопывание лент внутри каждого процесса!!!

    // !!!Цикл Reduc'ов в нулевой процесс
    std::vector<int> result_vector(number_of_matrix_rows);
    for (int i = 0; i < number_of_matrix_rows; ++i)
        MPI_Reduce(&result_vector_of_one_proc[i], &result_vector[i], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    // Цикл Reduc'ов в нулевой процесс!!!

    MPI_Type_free(&VerticalRibbonForOverloadProcsType);
    MPI_Type_free(&VerticalRibbonType);

    return result_vector;
}
