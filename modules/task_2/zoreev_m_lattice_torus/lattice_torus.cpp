// Copyright 2020 Zoreev Mikhail
#include "../../../modules/task_2/zoreev_m_lattice_torus/lattice_torus.h"

MPI_Comm createLatticeTorus(MPI_Comm comunicator, int dimensions_count) {
    int rank, size;
    int *dimensions = new int[dimensions_count]();
    int *periods = new int[dimensions_count];
    MPI_Comm result_communicator;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (int i = 0; i < dimensions_count; i++) {
        periods[i] = 1;
    }

    MPI_Dims_create(size, dimensions_count, dimensions);
    MPI_Cart_create(MPI_COMM_WORLD, dimensions_count, dimensions, periods, 0, &result_communicator);

    delete[] dimensions;
    delete[] periods;

    return result_communicator;
}

int LatticeTorusSend(void *buf, int count, MPI_Datatype datatype, int* dest_coords, int tag, MPI_Comm comm) {
    int dest;
    MPI_Cart_rank(comm, dest_coords, &dest);
    return MPI_Send(buf, count, datatype, dest, tag, comm);
}

int LatticeTorusRecv(void *buf, int count, MPI_Datatype datatype, int *source_coords, int tag, MPI_Comm comm,
                     MPI_Status *status) {
        int source;
        MPI_Cart_rank(comm, source_coords, &source);
        return MPI_Recv(buf, count, datatype, source, tag, comm, status);
    }
