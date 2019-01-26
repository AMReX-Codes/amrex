//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef MG_H
#define MG_H
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifndef MG_AGGLOMERATION_START
#define MG_AGGLOMERATION_START  8 // i.e. start the distributed v-cycle when boxes are smaller than 8^3
#endif
#ifndef MG_DEFAULT_BOTTOM_NORM
#define MG_DEFAULT_BOTTOM_NORM  1e-3
#endif
//------------------------------------------------------------------------------------------------------------------------------
typedef struct {
  int num_ranks;	// total number of MPI ranks for MPI_COMM_WORLD
  int my_rank;		// my MPI rank for MPI_COMM_WORLD
  int       num_levels;	// depth of the v-cycle
  level_type ** levels;	// array of pointers to levels

  struct {
    double MGBuild; // total time spent building the coefficients...
    double MGSolve; // total time spent in MGSolve
  }timers;
  int MGSolves_performed;
} mg_type;


//------------------------------------------------------------------------------------------------------------------------------
void          MGBuild(mg_type *all_grids, level_type *fine_grid, double a, double b, int minCoarseGridDim, const MPI_Comm comm);
void          MGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double dtol, double rtol);
void         FMGSolve(mg_type *all_grids, int onLevel, int u_id, int F_id, double a, double b, double rtol);
void            MGPCG(mg_type *all_grids, int onLevel, int x_id, int F_id, double a, double b, double dtol, double rtol);
void        MGDestroy(mg_type *all_grids);
void    MGPrintTiming(mg_type *all_grids, int fromLevel);
void    MGResetTimers(mg_type *all_grids);
void richardson_error(mg_type *all_grids, int levelh, int u_id);
//------------------------------------------------------------------------------------------------------------------------------
#endif
