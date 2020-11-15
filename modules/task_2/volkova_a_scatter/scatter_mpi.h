// Copyright 2020 Volkova Anastasia
#ifndef MODULES_TASK_2_VOLKOVA_A_SCATTER_SCATTER_MPI_H_
#define MODULES_TASK_2_VOLKOVA_A_SCATTER_SCATTER_MPI_H_


int MY_Scatter(void* sendBuf, int sendCount,
               MPI_Datatype sendType, void* recvBuf,
               int recvCount, MPI_Datatype recvType,
               int ROOT, MPI_Comm COMM);

#endif  // MODULES_TASK_2_VOLKOVA_A_SCATTER_SCATTER_MPI_H_
