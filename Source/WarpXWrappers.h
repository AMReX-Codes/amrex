#ifndef WARPX_WRAPPERS_H_
#define WARPX_WRAPPERS_H_

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

  void boxlib_init (int argc, char* argv[]);

  void boxlib_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);

  void boxlib_finalize (int finalize_mpi);

  void warpx_init ();

  void warpx_finalize ();

  void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.

#ifdef __cplusplus
}
#endif

#endif
