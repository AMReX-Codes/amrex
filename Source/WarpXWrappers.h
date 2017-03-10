#ifndef WARPX_WRAPPERS_H_
#define WARPX_WRAPPERS_H_

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

  void amrex_init (int argc, char* argv[]);

#ifdef BL_USE_MPI
  void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);
#endif

  void amrex_finalize (int finalize_mpi);

  void warpx_init ();

  void warpx_finalize ();

  void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.

  void addNParticles(int speciesnumber, int lenx, double* x, double* y, double* z, double* vx, double* vy, double* vz, int nattr, double* attr, int uniqueparticles);
  
  double warpx_getProbLo(int dir);

  double warpx_getProbHi(int dir);

  double* warpx_getArray();

#ifdef __cplusplus
}
#endif

#endif
