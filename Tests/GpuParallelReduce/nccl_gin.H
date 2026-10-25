#ifndef AMREX_TEST_NCCL_GIN_H_
#define AMREX_TEST_NCCL_GIN_H_

#include <nccl.h>
#include <mpi.h>

extern "C" {
void amrex_test_gin_initialize (int*, char***);
void amrex_test_gin_finalize ();
MPI_Comm amrex_test_gin_communicator ();
int amrex_test_gin_properties (ncclComm_t, int*, int*, int*);
int amrex_test_gin_open (ncclComm_t, int, void**);
int amrex_test_gin_run (void*, float const*, float*);
int amrex_test_gin_close (void*);
}

#endif
