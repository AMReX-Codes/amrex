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

  void addNParticles(int speciesnumber, int lenx, double* x, double* y, double* z, double* vx, double* vy, double* vz, int nattr, double* attr, int uniqueparticles);
  
  double warpx_getProbLo(int dir);

  double warpx_getProbHi(int dir);

#ifdef __cplusplus
}
#endif

#endif
