// Copyright 2020 Zoreev Mikhail
#ifndef MODULES_TASK_2_ZOREEV_M_LATTICE_TORUS_LATTICE_TORUS_H_
#define MODULES_TASK_2_ZOREEV_M_LATTICE_TORUS_LATTICE_TORUS_H_

#include <mpi.h>

MPI_Comm createLatticeTorus(MPI_Comm comunicator, int dimensions_count);

int LatticeTorusSend(void *buf, int count, MPI_Datatype datatype, int *dest_coords, int tag, MPI_Comm comm);

int LatticeTorusRecv(void *buf, int count, MPI_Datatype datatype, int *source_coords, int tag, MPI_Comm comm,
                     MPI_Status *status);

#endif  // MODULES_TASK_2_ZOREEV_M_LATTICE_TORUS_LATTICE_TORUS_H_
