// Copyright 2020 Napylov Evgenii

#include <iostream>
#include <vector>
#include <limits>
#include <string>
#include <random>
#include <ctime>
#include <mpi.h>
#include "../../../modules/task_2/napylov_e_gauss_horizontal/gauss_horizontal.h"

//#define DEBUG

const double EPSILON = std::numeric_limits<double>::epsilon();

// ------------------ DEBUG ------------------ //
void print_vec(std::vector<double> vec) {
    for (auto val : vec) {
        std::cout << val << ' ';
    }
    std::cout << std::endl;
}

void print_matrix(std::vector<double> vec, int rows, int cols) {
    for (auto i = 0; i < rows; i++) {
        for (auto j = 0; j < cols; j++) {
            std::cout << vec[i * cols + j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}
// ------------------ DEBUG ------------------ //

/*
    a1 a2 a3 | b1           
    a4 a5 a6 | b2   ->  a1 a2 a3 b1 a4 a5 a6 b2 a7 a8 a9 b3
    a7 a8 a9 | b3
*/

std::vector<double> RandomSystem(int rows, int cols) {
    srand(time(0));
    std::vector<double> result(rows * cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            result[i * cols + j] = ( ( double )rand() * ( 100 - 1 ) ) / ( double )RAND_MAX + 1;
        }
    }
    return result;
}

std::vector<double> SolveGaussSeq(std::vector<double> sys, int rows, int cols) {
    // ------------------ DEBUG ------------------ //
    #ifdef DEBUG
    print_matrix(sys, rows, cols);
    #endif
    // ------------------ DEBUG ------------------ //
    std::vector<double> result(rows);
    // прямой ход - исключение
    double exc_coeff = 0.0;
    for (auto curr_row = 0; curr_row < rows; curr_row++) {
        // на каждой интерации из нижних строк исключается переменная
        for (auto next_row = 1 + curr_row; next_row < rows; next_row++) {
            // расчитываем коэффициент, который исключит следующую переменную
            exc_coeff = sys[curr_row + cols * curr_row] / sys[curr_row + cols * next_row];
            // вычитаем из следующей строки * коэфф. текущую строку
            for (auto v_id = 0; v_id < cols; v_id++) {
                sys[v_id + cols * next_row] = exc_coeff * sys[v_id + cols * next_row] - sys[v_id + cols * curr_row];
            }
            // ------------------ DEBUG ------------------ //
            #ifdef DEBUG
            print_matrix(sys, rows, cols);
            #endif
            // ------------------ DEBUG ------------------ //
        }
    }
    // ------------------ DEBUG ------------------ //
    #ifdef DEBUG
    print_matrix(sys, rows, cols);
    #endif
    // ------------------ DEBUG ------------------ //

    // обратный ход - вычисление
    // проходим снизу вверх и вычисляем неизвестные, подставляя уже найденные
    for (auto curr_row = rows - 1; curr_row > -1; curr_row--) {
        // это правая часть системы
        result[curr_row] = sys[curr_row * cols + cols - 1];
        // вычитаем из правой части уже найденные * на их коэфф.
        for (auto prev_row = rows - 1; prev_row > curr_row; prev_row--) {
            result[curr_row] -= sys[curr_row * cols + prev_row] * result[prev_row];
        }
        // находим неизв. (избавляемся от коэфф.)
        result[curr_row] /= sys[curr_row * cols + curr_row];
    }
    return result;
}

bool CheckSolution(std::vector<double> sys, int rows, int cols, std::vector<double> answer, double epsilon) {
    double tmp_sum = 0.0;
    int id;
    int id_ans;
    for (int row = 0; row < rows; row++) {
        id_ans = 0;
        for (int col = 0; col < cols - 1; col++) {
            tmp_sum += (sys[row * cols + col] * answer[id_ans]);
            id = row * cols + col;
            id_ans++;
        }
        if ((tmp_sum - sys[id + 1] >= epsilon) || (tmp_sum - sys[id + 1] <= -epsilon)) {
            return false;
        }
        tmp_sum = 0.0;
    }
    return true;
}

/*std::vector<double>*/ void SolveGaussParallel(std::vector<double> sys, int rows, int cols) {
    // -------------------- Подготовка -------------------- //
    //Создаем новый коммуникатор, который исключает лишние процессы (world_size > rows)
    MPI_Group group_world;
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int work_proc_size;
    size > rows ? work_proc_size = rows : work_proc_size = size;
    int *work_rank_array = new int[work_proc_size];
    for (int i = 0; i < work_proc_size; i++) work_rank_array[i] = i;
    MPI_Group group_work;
    MPI_Comm COMM_WORK;
    MPI_Group_incl(group_world, work_proc_size, work_rank_array, &group_work);
    MPI_Comm_create(MPI_COMM_WORLD, group_work, &COMM_WORK);
    // -------------------- Конец подготовки -------------------- //

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    MPI_Barrier(MPI_COMM_WORLD);
    if (world_rank < work_proc_size) {
        int rank;
        const int root = 0;
        MPI_Comm_rank(COMM_WORK, &rank);
        MPI_Comm_size(COMM_WORK, &work_proc_size);
        const unsigned int rows_in_process = rows / size;
        const unsigned int rows_rem = rows % size;

        // раскладываем поровну + остатки "вторым проходом"
        const unsigned int local_vec_size = (rows_in_process + (rank < rows_rem ? 1 : 0)) * cols;

        std::vector<double> local_vec(local_vec_size);

        if (rank == 0) {
            for (int proc = 1; proc < work_proc_size; proc++) {
                int tmp_size = (rows_in_process + (proc < rows_rem ? 1 : 0)) * cols;
                MPI_Send(sys.data() + proc * tmp_size, tmp_size, MPI_DOUBLE, proc, 0, COMM_WORK);
            }
            local_vec = std::vector<double>(sys.begin(), sys.begin() + local_vec_size);
        } else {
            MPI_Status status;
            MPI_Recv(local_vec.data(), local_vec_size, MPI_DOUBLE, 0, 0, COMM_WORK, &status);
        }
        // ждем, пока данные появятся во всех процессах
        // пока не понятно насколько это необходимо
        MPI_Barrier(COMM_WORK);

        // ТЕПЕРЬ ДАННЫЕ РАСПРЕДЕЛЕНЫ

        std::cout << rank << ": ";
        print_vec(local_vec);
    }
}
