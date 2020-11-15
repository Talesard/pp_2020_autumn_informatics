// Copyright 2020 Gogov Vladislav
#include <mpi.h>
#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <limits>
#include "../../../modules/task_2/gogov_v_block_striped_columnwise/block_striped_columnwise.h"

bool assertVectors(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size())
        throw "Different size";
    for (size_t i = 0; i < vec1.size(); i++) {
        if ((std::fabs(vec1[i] - vec2[i]) >= std::numeric_limits<double>::epsilon() * 1000000.0))
            return false;
    }
    return true;
}

void printVector(const std::vector<double>& vec, int vecSize) {
    std::cout << "Vector print:" << std::endl;
    for (int i = 0; i < vecSize; i++)
        std::cout << vec[i] << std::endl;
}

void printMatrix(const std::vector<double>& vec, int rows, int cols) {
    std::cout << "Matrix print:" << std::endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            std::cout << vec[i * cols + j] << "  ";
        }
        std::cout << std::endl;
    }
}

Matrix createRandomMatrix(int rows, int cols) {
    std::random_device rd;
    std::mt19937 mersenne(rd());
    std::uniform_real_distribution<> urd(-50.0 , 50.0);
    Matrix matrix(rows * cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            matrix[i * cols + j] = static_cast<double>(urd(mersenne));
        }
    }
    return matrix;
}

std::vector<double> multMatrixByVectorSequential(const Matrix& matrix, int rows, int cols,
                                            const std::vector<double>& vec, int vecSize) {
    if (cols != vecSize) {
        throw("Invalid size");
    }
    std::vector<double> result(rows, 0);
    for (int i = 0; i < cols; i++) {
        for (int j = 0; j < rows; j++) {
            result[j] += matrix[j * cols + i] * vec[i];
        }
    }
    return result;
}

void arraysDistribution(int* sendNum, int* sendInd, int n, int procNum, int remain, int cols) {
    sendNum[0] = (remain != 0 ? n + 1 : n);
    sendInd[0] = 0;
    for (int i = 1; i < procNum; i++) {
        sendNum[i] = (i < remain ? n + 1 : n);
        sendInd[i] = sendInd[i - 1] + sendNum[i - 1];
    }
}

std::vector<double> sumMatrixByCols(const std::vector<double>& vec, int rows, int cols) {
    std::vector<double> result(cols, 0);
    for (int i = 0; i < cols; i++) {
        double temp = 0;
        for (int j = 0; j < rows; j++) {
            temp += vec[j * cols + i];
        }
        result[i] = temp;
    }
    return result;
}

std::vector<double> multMatrixByVectorParallel(const Matrix& matrix, int rows, int cols,
                                            const std::vector<double>& vec, int vecSize) {
    if (cols != vecSize) {
        throw("Invalid size");
    }
    int procNum, procRank;
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
    int n = cols / procNum;
    int remain = cols - n * procNum;

    //  Определение массивов для type
    std::vector<int> sendNum(procNum);
    std::vector<int> sendInd(procNum);
    arraysDistribution(sendNum.data(), sendInd.data(), n, procNum, remain, cols);

    //  Определение производного типа столбец
    int countType = (remain == 0 ? 1 : 2);
    MPI_Datatype* colType = new MPI_Datatype[countType];
    for (int i = 0; i < countType; i++) {
        MPI_Type_vector(rows, sendNum[i == 0 ? i : procNum - 1], cols, MPI_DOUBLE, &colType[i]);
        MPI_Type_commit(&colType[i]);
    }

    //  Отправка вектора
    std::vector<double> partVector(sendNum[procRank]);
    MPI_Scatterv(vec.data(), sendNum.data(), sendInd.data(), MPI_DOUBLE,
                partVector.data(), sendNum[procRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //  Передача и прием матрицы
    std::vector<double> colVector(sendNum[procRank] * rows);
    if (procRank == 0) {
        for (int i = 1; i < procNum; i++) {
            MPI_Send(matrix.data() + sendInd[i], 1, colType[remain == 0 || i < remain ? 0 : 1], i, 0, MPI_COMM_WORLD);
        }
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < sendNum[procRank]; j++) {
                colVector[i * sendNum[procRank] + j] = matrix[i * cols + j];
            }
        }
    } else {
        MPI_Status status;
        MPI_Recv(colVector.data(), sendNum[procRank] * rows, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    }
    std::vector<double> localVector = multMatrixByVectorSequential(colVector, rows, sendNum[procRank],
                                                                partVector, sendNum[procRank]);

    //  Сбор промежуточных векторов
    std::vector<double> tempMatrix(procNum * rows);
    MPI_Gather(localVector.data(), static_cast<int>(localVector.size()), MPI_DOUBLE,
                tempMatrix.data(), rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Суммирование столбцов
    std::vector<double> resultVector(0);
    if (procRank == 0) {
        resultVector = sumMatrixByCols(tempMatrix, procNum, rows);
    }

    for (int i = 0; i < countType; i++)
        MPI_Type_free(&colType[i]);
    return resultVector;
}
