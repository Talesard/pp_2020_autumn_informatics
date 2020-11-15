// Copyright 2020 Sidorova Alexandra
#ifndef MODULES_TASK_2_SIDOROVA_A_ALLREDUCE_ALLREDUCE_H_
#define MODULES_TASK_2_SIDOROVA_A_ALLREDUCE_ALLREDUCE_H_

#include <mpi.h>

int Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm);

#endif  // MODULES_TASK_2_SIDOROVA_A_ALLREDUCE_ALLREDUCE_H_
