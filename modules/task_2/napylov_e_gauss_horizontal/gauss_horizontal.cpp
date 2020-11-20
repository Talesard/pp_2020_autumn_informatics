// Copyright 2020 Napylov Evgenii
#include "../../../modules/task_2/napylov_e_gauss_horizontal/gauss_horizontal.h"
#include <mpi.h>
#include <iostream>
#include <vector>

// #define DEBUG_PAR
// #define DEBUG_SEQ


void print_vec(std::vector<double> vec) {
    for (auto val : vec) {
        std::cout << val << ' ';
    }
    std::cout << std::endl;
}

void print_vec(std::vector<int> vec) {
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

/*
    a1 a2 a3 | b1
    a4 a5 a6 | b2   ->  a1 a2 a3 b1 a4 a5 a6 b2 a7 a8 a9 b3
    a7 a8 a9 | b3
*/

std::vector<double> SystemForPerformanceTest(int rows, int cols) {
    std::vector<double> res(rows * cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            res[i * cols + j] = 5;
        }
    }
    return res;
}

std::vector<double> SolveGaussSeq(std::vector<double> sys, int rows, int cols) {
    #ifdef DEBUG_SEQ
    print_matrix(sys, rows, cols);
    #endif
    std::vector<double> result(rows);
    // прямой ход - исключение
    double exc_coeff = 0.0;
    for (int curr_r = 0; curr_r < rows; curr_r++) {
        // на каждой интерации из нижних строк исключается переменная
        for (int next_r = 1 + curr_r; next_r < rows; next_r++) {
            // расчитываем коэффициент, который исключит следующую переменную
            exc_coeff = sys[curr_r + cols * curr_r] / sys[curr_r + cols * next_r];
            // вычитаем из следующей строки * коэфф. текущую строку
            for (int v_id = 0; v_id < cols; v_id++) {
                sys[v_id + cols * next_r] = exc_coeff * sys[v_id + cols * next_r] - sys[v_id + cols * curr_r];
            }
            #ifdef DEBUG_SEQ
            print_matrix(sys, rows, cols);
            #endif
        }
    }
    #ifdef DEBUG_SEQ
    print_matrix(sys, rows, cols);
    #endif

    // обратный ход - вычисление
    // проходим снизу вверх и вычисляем неизвестные, подставляя уже найденные
    for (int curr_r = rows - 1; curr_r > -1; curr_r--) {
        // вычитаем уже найденные * на их коэфф.
        for (int prev_r = rows - 1; prev_r > curr_r; prev_r--) {
            result[curr_r] -= sys[curr_r * cols + prev_r] * result[prev_r];
        }
        // + свободный член
        result[curr_r] += sys[curr_r * cols + cols - 1];
        // находим неизв. (избавляемся от коэфф.)
        result[curr_r] /= sys[curr_r * cols + curr_r];
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

std::vector<double> SolveGaussParallel(std::vector<double> sys, int rows, int cols) {
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
    MPI_Barrier(MPI_COMM_WORLD);
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::vector<double> result(rows);

    if (world_rank < work_proc_size) {
        int rank;
        MPI_Comm_rank(COMM_WORK, &rank);
        MPI_Comm_size(COMM_WORK, &work_proc_size);
        const int rows_in_process = rows / size;
        const int rows_rem = rows % size;
        // раскладываем поровну + 1 из остатка
        int local_vec_size;
        if (rank < rows_rem) {
            local_vec_size = (rows_in_process + 1) * cols;
        } else {
            local_vec_size = rows_in_process * cols;
        }
        std::vector<double> local_vec(local_vec_size);
        std::vector<int> map(rows);  // map[row]=proc, в котором лежит row
        std::vector<int> size_vec;  // размеры для gatherv
        int *displ;  // смещения для gatherv

        if (rank == 0) {
            displ = new int[size];
            displ[0] = 0;
            size_vec.resize(size);
            local_vec = std::vector<double>(sys.begin(), sys.begin() + local_vec_size);
            int offset = local_vec_size;
            for (int i = 0; i < local_vec_size / cols; i++) {
                map[i] = 0;
            }
            size_vec[0] = local_vec_size;
            for (int proc = 1; proc < work_proc_size; proc++) {
                int tmp_size;
                if (proc < rows_rem) {
                    tmp_size = (rows_in_process + 1) * cols;
                } else {
                    tmp_size = rows_in_process * cols;
                }
                size_vec[proc] = tmp_size;
                displ[proc] = tmp_size;
                MPI_Send(sys.data() + offset, tmp_size, MPI_DOUBLE, proc, 0, COMM_WORK);
                for (int i = offset / cols; i < offset / cols + tmp_size / cols; i++) {
                    map[i] = proc;
                }
                offset += tmp_size;
            }
            for (int i = 1; i < size; i++) displ[i] = displ[i-1] + size_vec[i-1];

            #ifdef DEBUG_PAR
            std::cout << "map: "; print_vec(map);
            #endif
        } else {
            MPI_Status status;
            MPI_Recv(local_vec.data(), local_vec_size, MPI_DOUBLE, 0, 0, COMM_WORK, &status);
        }
        MPI_Bcast(map.data(), rows, MPI_INT, 0, COMM_WORK);

        int local_rows_count = local_vec_size / cols;
        /*
            local_vec - локальный вектор, содержащий все строки процесса
            local_vec_size - его размер
            local_rows_count - число строк в процессе
            map - в каком процессе строка row: map[row] = proc
        */
        int local_id = 0;
        // рабочая строка, с помощью которой процессы будут исключать
        std::vector<double> stroka(cols);
        for (int r = 0; r < rows - 1; r++) {
            int proc = map[r];
            if (rank == proc) {
                stroka = std::vector<double>(local_vec.data() + local_id * cols,
                                            local_vec.data() + local_id * cols + cols);
                local_id++;
                for (int pn = rank + 1; pn < work_proc_size; pn++) {
                    MPI_Send(stroka.data(), cols, MPI_DOUBLE, pn, 0, MPI_COMM_WORLD);
                }
            }
            if (rank > map[r]) {
                MPI_Status st;
                MPI_Recv(stroka.data(), cols, MPI_DOUBLE, map[r], 0, MPI_COMM_WORLD, &st);
            }

            // прямой ход - исключение
            if (rank == proc) {
                // local_id не уведичиваем на 1, т.к. это сделано выше
                for (int l_r = local_id; l_r < local_rows_count; l_r++) {
                    double coeff = stroka[r] / local_vec[l_r * cols + r];
                    for (int el_id = 0; el_id < cols; el_id++) {
                        local_vec[l_r * cols + el_id] = local_vec[l_r * cols + el_id] * coeff - stroka[el_id];
                    }
                }
            } else {
                if (rank > proc) {
                    for (int l_r = 0; l_r < local_rows_count; l_r++) {
                        double coeff = stroka[r] / local_vec[l_r * cols + r];
                        for (int el_id = 0; el_id < cols; el_id++) {
                            local_vec[l_r * cols + el_id] = local_vec[l_r * cols + el_id] * coeff - stroka[el_id];
                        }
                    }
                }
            }
        }
        #ifdef DEBUG_PAR
        std::cout << "rank: " << rank << ", size: " << local_vec_size << ", vec: "; print_vec(local_vec);
        #endif

        // ТЕПЕРЬ ВСЕ СТРОКИ ВСЕХ ПРОЦЕССОВ СОСТАВЛЯЮТ ВЕРХНЕ-ТРЕУГОЛЬНУЮ МАТРИЦУ

        // обратный ход - исключение
        // т.к. вычисление обратным ходом в любом случае последовательное, пусть этим занимается нулевой прцесс
        for (int p = 1; p < size; p++) {
            if (rank == p) {
                MPI_Send(local_vec.data(), local_vec_size, MPI_DOUBLE, 0, p, MPI_COMM_WORLD);
            }
        }

        if (rank == 0) {
            MPI_Status st;
            for (int i = 0; i < local_vec_size; i++) {
                sys[i] = local_vec[i];
            }
            for (int p = 1; p < work_proc_size; p++) {
                MPI_Recv(sys.data() + displ[p], size_vec[p], MPI_DOUBLE, p, p, MPI_COMM_WORLD, &st);
            }
        }
        MPI_Barrier(COMM_WORK);
        if (rank == 0) {
            #ifdef DEBUG_PAR
            print_matrix(sys, rows, cols);
            #endif
            for (int curr_r = rows - 1; curr_r > -1; curr_r--) {
                result[curr_r] = 0.0;
                for (int prev_r = rows - 1; prev_r > curr_r; prev_r--) {
                    result[curr_r] -= sys[curr_r * cols + prev_r] * result[prev_r];
                }
                result[curr_r] += sys[curr_r * cols + cols - 1];
                result[curr_r] /= sys[curr_r * cols + curr_r];
            }
        }
        MPI_Comm_free(&COMM_WORK);
    }
    // можно не отправлять всем, но тогда правильный ответ вернется только в нулевом процессе
    // MPI_Bcast(result.data(), rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Group_free(&group_work);
    return result;
}
