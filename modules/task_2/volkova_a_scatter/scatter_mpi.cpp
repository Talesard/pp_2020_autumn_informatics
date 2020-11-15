// Copyright 2020 Volkova Anastasia
#include <mpi.h>
#include <cstring>
#include <iostream>
#include "../../../modules/task_2/volkova_a_scatter/scatter_mpi.h"

int MY_Scatter(void* sendBuf, int sendCount,
    MPI_Datatype sendType, void* recvBuf,
    int recvCount, MPI_Datatype recvType,
    int ROOT, MPI_Comm COMM) {
    if (sendCount != recvCount)return MPI_ERR_COUNT;
    if (sendCount <= 0) return MPI_ERR_COUNT;
    if (recvCount <= 0) return MPI_ERR_COUNT;
    if (sendBuf == nullptr) return MPI_ERR_BUFFER;
    int size;
    int rank;
    int sizeTypeSend;
    int tag = 15;
    if (COMM == MPI_COMM_NULL) {
        return MPI_ERR_COMM;
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        MPI_Status stat;
        return stat.MPI_ERROR;
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        MPI_Status stat;
        return stat.MPI_ERROR;
    }
    if (MPI_Type_size(sendType, &sizeTypeSend) != MPI_SUCCESS) {
        MPI_Status stat;
        return stat.MPI_ERROR;
    }

    if (ROOT >= size || ROOT < 0)return MPI_ERR_ROOT;

    int sizeTypeCount = sendCount * sizeTypeSend;
    if (rank == ROOT) {
        for (int proc = 0; proc < size; ++proc) {
            if (proc != ROOT) {
                int pos = proc * sizeTypeCount;
                MPI_Send(static_cast<char*>(sendBuf) + pos, sendCount, sendType, proc, tag, COMM);
            }
        }
    }
    if (rank == ROOT) {
        int k = rank * sizeTypeCount;
        std::memcpy(recvBuf, static_cast<char*>(sendBuf) + k, sizeTypeCount);
    } else {
        MPI_Status STATUS;
        MPI_Recv(recvBuf, recvCount, recvType, ROOT, tag, COMM, &STATUS);
    }
    return MPI_SUCCESS;
}


