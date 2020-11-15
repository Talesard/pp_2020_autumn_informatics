// Copyright 2020 Sidorova Alexandra

#include <mpi.h>
#include <memory>
#include <cstring>
#include "../../../modules/task_2/sidorova_a_allreduce/allreduce.h"

int Allreduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {
    if (sendbuf == nullptr || recvbuf == nullptr)
        return MPI_ERR_BUFFER;
    if (count <= 0)
        return MPI_ERR_COUNT;
    if (op != MPI_MAX && op != MPI_MIN && op != MPI_SUM && op != MPI_PROD)
        return MPI_ERR_OP;
    if (comm != MPI_COMM_WORLD && comm != MPI_COMM_SELF)
        return MPI_ERR_COMM;

    int size, rank;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    void* buffer;

    if (datatype == MPI_INT) {
        std::memcpy(recvbuf, sendbuf, count * sizeof(int));
        buffer = new int[count];
    } else if (datatype == MPI_FLOAT) {
        std::memcpy(recvbuf, sendbuf, count * sizeof(float));
        buffer = new float[count];
    } else if (datatype == MPI_DOUBLE) {
        std::memcpy(recvbuf, sendbuf, count * sizeof(double));
        buffer = new double[count];
    } else {
        return MPI_ERR_TYPE;
    }

    int localSize = size;

    while (localSize != 1) {
        int half = localSize / 2;

        for (int i = 0; i < half; ++i)
            if (rank == localSize - i - 1)
                MPI_Send(recvbuf, count, datatype, i, 0, comm);
        for (int i = 0; i < half; ++i) {
            if (rank == i) {
                MPI_Status status;
                if (datatype == MPI_INT) {
                    MPI_Recv(buffer, count, datatype, localSize - i - 1, 0, comm, &status);

                    if (op == MPI_MAX) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<int *>(recvbuf)[j] < static_cast<int *>(buffer)[j])
                                static_cast<int *>(recvbuf)[j] = static_cast<int *>(buffer)[j];
                        }
                    } else if (op == MPI_MIN) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<int *>(recvbuf)[j] > static_cast<int *>(buffer)[j])
                                static_cast<int *>(recvbuf)[j] = static_cast<int *>(buffer)[j];
                        }
                    } else if (op == MPI_SUM) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<int*>(recvbuf)[j] += static_cast<int*>(buffer)[j];
                        }
                    } else if (op == MPI_PROD) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<int *>(recvbuf)[j] *= static_cast<int *>(buffer)[j];
                        }
                    }
                } else if (datatype == MPI_FLOAT) {
                    MPI_Recv(buffer, count, datatype, localSize - i - 1, 0, comm, &status);

                    if (op == MPI_MAX) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<float *>(recvbuf)[j] < static_cast<float *>(buffer)[j])
                                static_cast<float *>(recvbuf)[j] = static_cast<float *>(buffer)[j];
                        }
                    } else if (op == MPI_MIN) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<float *>(recvbuf)[j] > static_cast<float *>(buffer)[j])
                                static_cast<float *>(recvbuf)[j] = static_cast<float *>(buffer)[j];
                        }
                    } else if (op == MPI_SUM) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<float *>(recvbuf)[j] += static_cast<float *>(buffer)[j];
                        }
                    } else if (op == MPI_PROD) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<float *>(recvbuf)[j] *= static_cast<float *>(buffer)[j];
                        }
                    }
                } else if (datatype == MPI_DOUBLE) {
                    MPI_Recv(buffer, count, datatype, localSize - i - 1, 0, comm, &status);

                    if (op == MPI_MAX) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<double *>(recvbuf)[j] < static_cast<double *>(buffer)[j])
                                static_cast<double *>(recvbuf)[j] = static_cast<double *>(buffer)[j];
                        }
                    } else if (op == MPI_MIN) {
                        for (int j = 0; j < count; ++j) {
                            if (static_cast<double *>(recvbuf)[j] > static_cast<double *>(buffer)[j])
                                static_cast<double *>(recvbuf)[j] = static_cast<double *>(buffer)[j];
                        }
                    } else if (op == MPI_SUM) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<double *>(recvbuf)[j] += static_cast<double *>(buffer)[j];
                        }
                    } else if (op == MPI_PROD) {
                        for (int j = 0; j < count; ++j) {
                            static_cast<double *>(recvbuf)[j] *= static_cast<double *>(buffer)[j];
                        }
                    }
                }
            }
        }

        localSize = half + localSize % 2;
    }

    if (size != 1) {
        if (rank == 0) {
            for (int i = 1; i < size; ++i) {
                MPI_Send(recvbuf, count, datatype, i, i, comm);
            }
        } else {
            MPI_Status status;
            MPI_Recv(recvbuf, count, datatype, 0, rank, comm, &status);
        }
    }

    if (datatype == MPI_INT)
        delete[] static_cast<int*>(buffer);
    else if (datatype == MPI_FLOAT)
        delete[] static_cast<float*>(buffer);
    else if (datatype == MPI_DOUBLE)
        delete[] static_cast<double*>(buffer);

    return MPI_SUCCESS;
}
