// Copyright 2020 Napylov Evgenii

#include <iostream>
#include <vector>
#include <limits>
#include <string>
#include <mpi.h>
#include "../../../modules/task_2/napylov_e_gauss_horizontal/gauss_horizontal.h"

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

std::vector<double> SolveGaussSeq(std::vector<double> sys, int rows, int cols) {
    // ------------------ DEBUG ------------------ //
    bool debug = false;
    if (debug) print_matrix(sys, rows, cols);
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
            if (debug) print_matrix(sys, rows, cols);
            // ------------------ DEBUG ------------------ //
        }
    }
    // ------------------ DEBUG ------------------ //
    if (debug) print_matrix(sys, rows, cols);
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
    bool debug = false;
    double tmp_sum = 0.0;
    int id;
    int id_ans;
    for (int row = 0; row < rows; row++) {
        id_ans = 0;
        for (int col = 0; col < cols - 1; col++) {
            if (debug) std::cout << sys[row * cols + col] << " * " << answer[id_ans] << std::endl;
            tmp_sum += (sys[row * cols + col] * answer[id_ans]);
            id = row * cols + col;
            id_ans++;
        }
        if (debug) std::cout << "tmp_sum = " << tmp_sum - sys[id + 1] << std::endl;
        if ((tmp_sum - sys[id + 1] >= epsilon) || (tmp_sum - sys[id + 1] <= -epsilon)) {
            if (debug) std::cout << "!!!" << std::endl;
            return false;
        }
        tmp_sum = 0.0;
        if (debug) std::cout << std::endl;
    }
    return true;
}

std::vector<double> SolveGaussParallel(std::vector<double> sys, int rows, int cols) {
    bool debug = false;
    // sys уже есть во всех процессах
    int size, rank;
    const int root = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // тут нужно будет создать новый комм если процессов > строк.   upd: решилось САМО Pog
    // и подумать что делать если <
    // ...

    // Отправляем процессу i строку i
    std::vector<double> local_str(cols);
    MPI_Scatter(sys.data(), cols, MPI_DOUBLE, local_str.data(), cols, MPI_DOUBLE, root, MPI_COMM_WORLD);

    // прямой ход
    std::vector<double> stroka(cols);
    for (int i = 0; i < rows - 1; i++) {
        if (rank == i) {
            stroka = local_str;
        }
        // процесс i передает свою строку i всем
        MPI_Bcast(stroka.data(), cols, MPI_DOUBLE, i, MPI_COMM_WORLD);

        // ------------------ DEBUG ------------------ //
        if (debug) {
            MPI_Barrier(MPI_COMM_WORLD); 
            std::cout << "i=" << i << " r="<< rank << ": ";
            print_vec(stroka);
            std::cout << "--------" << std::endl;
            MPI_Barrier(MPI_COMM_WORLD);
        }
        // ------------------ DEBUG ------------------ //

        // процессы i+1...rows исключают из совей строки Xi
        for (int j = i + 1; j < rows; j++) {
            if (rank == j) {
                // тут исключение Xi из local_str используя принятую из bcast stroka 
                double coeff = stroka[i] / local_str[i];
                // local_str*coeff - stroka
                for (int k = 0; k < cols; k++) {
                    local_str[k] = local_str[k] * coeff - stroka[k];
                }
                // ------------------ DEBUG ------------------ //
                if (debug) print_vec(local_str);
                if (debug) std::cout << "-------" << std::endl;
                // ------------------ DEBUG ------------------ //
            }
        } // кажется, это реально работает... Pog?
    } // конец прямого хода
    // Теперь local_str всех процессов образуют верхне-треугольную матрицу



    // обратный ход - вычисление неизвестных
    // Это работает последовательно KEKWait. 
    // Подумать, как можно распараллелить
    std::vector<double> result(rows);
    for (int i = rows - 1; i > -1; i--) { // цикл по процессам
        result[i] = local_str[cols - 1]; // свободный член, из него будем вычитать все, что слева
        for (int j = rows - 1; j > i; j--) {
            result[i] -= local_str[j] * result[j];
        }
        result[i] /= local_str[i];
        // отправляем результат всем процессам. (а может нужно отправить только процессу #i-1 ???)
        // и отправлять только существующие значения, а не весь вектор, где много 0
        MPI_Bcast(result.data(), rows, MPI_DOUBLE, i, MPI_COMM_WORLD);
    }

    return result;
}
