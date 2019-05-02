#ifndef WARPX_WRAPPERS_H_
#define WARPX_WRAPPERS_H_

#ifdef BL_USE_MPI
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

    int warpx_nSpecies();

    bool warpx_use_fdtd_nci_corr();

    int warpx_l_lower_order_in_v();

    int warpx_nComps();

    int warpx_SpaceDim();

    void amrex_init (int argc, char* argv[]);

#ifdef BL_USE_MPI
    void amrex_init_with_inited_mpi (int argc, char* argv[], MPI_Comm mpicomm);
#endif

    void amrex_finalize (int finalize_mpi);
    
    void warpx_init ();
    
    void warpx_finalize ();

    typedef void(*WARPX_CALLBACK_PY_FUNC_0)();

    void warpx_set_callback_py_afterinit (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforeEsolve (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterEsolve (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforedeposition (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterdeposition (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particlescraper (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particleloader (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_beforestep (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterstep (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_afterrestart (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_particleinjection (WARPX_CALLBACK_PY_FUNC_0);
    void warpx_set_callback_py_appliedfields (WARPX_CALLBACK_PY_FUNC_0);

    void warpx_evolve (int numsteps);  // -1 means the inputs parameter will be used.
    
    void warpx_addNParticles(int speciesnumber, int lenx,
                             double* x, double* y, double* z,
                             double* vx, double* vy, double* vz,
                             int nattr, double* attr, int uniqueparticles);

    void warpx_ConvertLabParamsToBoost();
  
    double warpx_getProbLo(int dir);
    
    double warpx_getProbHi(int dir);
    
    long warpx_getNumParticles(int speciesnumber);
    
    double** warpx_getEfield(int lev, int direction, 
                             int *return_size, int* ngrow, int **shapes);
    
    int* warpx_getEfieldLoVects(int lev, int direction, 
                                int *return_size, int* ngrow);
    
    double** warpx_getBfield(int lev, int direction, 
                             int *return_size, int* ngrow, int **shapes);
    
    int* warpx_getBfieldLoVects(int lev, int direction, 
                                int *return_size, int* ngrow);
    
    double** warpx_getCurrentDensity(int lev, int direction, 
                                     int *return_size, int* ngrow, int **shapes);
    
    int* warpx_getCurrentDensityLoVects(int lev, int direction, 
                                        int *return_size, int* ngrow);
    
    double** warpx_getParticleStructs(int speciesnumber,
                                      int* num_tiles, int** particles_per_tile);
    
    double** warpx_getParticleArrays(int speciesnumber, int comp,
                                     int* num_tiles, int** particles_per_tile);

  void warpx_ComputeDt ();
  void warpx_MoveWindow ();

  void warpx_EvolveE (double dt);
  void warpx_EvolveB (double dt);
  void warpx_FillBoundaryE ();
  void warpx_FillBoundaryB ();
  void warpx_SyncCurrent ();
  void warpx_UpdateAuxilaryData ();
  void warpx_PushParticlesandDepose (double cur_time);

  int warpx_getistep (int lev);
  void warpx_setistep (int lev, int ii);
  double warpx_gett_new (int lev);
  void warpx_sett_new (int lev, double time);
  double warpx_getdt (int lev);

  int warpx_maxStep ();
  double warpx_stopTime ();

  int warpx_checkInt ();
  int warpx_plotInt ();

  void warpx_WriteCheckPointFile ();
  void warpx_WritePlotFile ();

  int warpx_finestLevel ();

  void mypc_Redistribute ();

#ifdef __cplusplus
}
#endif

#endif
