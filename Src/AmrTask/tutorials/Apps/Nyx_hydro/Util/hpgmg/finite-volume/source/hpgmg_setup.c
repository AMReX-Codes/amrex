//------------------------------------------------------------------------------------------------------------------------------
// Copyright Notice
//------------------------------------------------------------------------------------------------------------------------------
// HPGMG, Copyright (c) 2014, The Regents of the University of
// California, through Lawrence Berkeley National Laboratory (subject to
// receipt of any required approvals from the U.S. Dept. of Energy).  All
// rights reserved.
//
// If you have questions about your rights to use or distribute this
// software, please contact Berkeley Lab's Technology Transfer Department
// at  TTD@lbl.gov.
//
// NOTICE.  This software is owned by the U.S. Department of Energy.  As
// such, the U.S. Government has been granted for itself and others
// acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
// license in the Software to reproduce, prepare derivative works, and
// perform publicly and display publicly.  Beginning five (5) years after
// the date permission to assert copyright is obtained from the U.S.
// Department of Energy, and subject to any subsequent five (5) year
// renewals, the U.S. Government is granted for itself and others acting
// on its behalf a paid-up, nonexclusive, irrevocable, worldwide license
// in the Software to reproduce, prepare derivative works, distribute
// copies to the public, perform publicly and display publicly, and to
// permit others to do so.
//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
//------------------------------------------------------------------------------------------------------------------------------
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
//------------------------------------------------------------------------------------------------------------------------------
#include "defines.h"
#include "level.h"
#include "mg.h"
#include "operators.h"
#include "solvers.h"
//------------------------------------------------------------------------------------------------------------------------------
void hpgmg_setup(const int log2_box_dim,
                 const int target_boxes_per_rank,
                 const int OMP_Threads,
                 const int OMP_Nested,
                 const int requested_threading_model,
                 const int actual_threading_model) {
  int my_rank=0;
  int num_tasks=1;

  #ifdef USE_MPI
  MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//if(actual_threading_model>requested_threading_model)actual_threading_model=requested_threading_model;
  if(my_rank==0){
       if(requested_threading_model == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else if(requested_threading_model == MPI_THREAD_SINGLE    )printf("Requested MPI_THREAD_SINGLE, ");
  else if(requested_threading_model == MPI_THREAD_FUNNELED  )printf("Requested MPI_THREAD_FUNNELED, ");
  else if(requested_threading_model == MPI_THREAD_SERIALIZED)printf("Requested MPI_THREAD_SERIALIZED, ");
  else if(requested_threading_model == MPI_THREAD_MULTIPLE  )printf("Requested MPI_THREAD_MULTIPLE, ");
  else                                                       printf("Requested Unknown MPI Threading Model (%d), ",requested_threading_model);
       if(actual_threading_model    == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else if(actual_threading_model    == MPI_THREAD_SINGLE    )printf("got MPI_THREAD_SINGLE\n");
  else if(actual_threading_model    == MPI_THREAD_FUNNELED  )printf("got MPI_THREAD_FUNNELED\n");
  else if(actual_threading_model    == MPI_THREAD_SERIALIZED)printf("got MPI_THREAD_SERIALIZED\n");
  else if(actual_threading_model    == MPI_THREAD_MULTIPLE  )printf("got MPI_THREAD_MULTIPLE\n");
  else                                                       printf("got Unknown MPI Threading Model (%d)\n",actual_threading_model);
  }
  #endif


  if(log2_box_dim<4){
    if(my_rank==0){printf("log2_box_dim must be at least 4\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(target_boxes_per_rank<1){
    if(my_rank==0){printf("target_boxes_per_rank must be at least 1\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }

  if(my_rank==0){
    if(OMP_Nested)fprintf(stdout,"%d MPI Tasks of %d threads (OMP_NESTED=TRUE)\n\n" ,num_tasks,OMP_Threads);
             else fprintf(stdout,"%d MPI Tasks of %d threads (OMP_NESTED=FALSE)\n\n",num_tasks,OMP_Threads);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // calculate the problem size...
  #ifndef MAX_COARSE_DIM
  #define MAX_COARSE_DIM 11
  #endif
  int64_t box_dim=1<<log2_box_dim;
  int64_t target_boxes = (int64_t)target_boxes_per_rank*(int64_t)num_tasks;
  int64_t boxes_in_i = -1;
  int64_t bi;
  for(bi=1;bi<1000;bi++){ // all possible problem sizes
    int64_t total_boxes = bi*bi*bi;
    if(total_boxes<=target_boxes){
      int64_t coarse_grid_dim = box_dim*bi;
      while( (coarse_grid_dim%2) == 0){coarse_grid_dim=coarse_grid_dim/2;}
      if(coarse_grid_dim<=MAX_COARSE_DIM){
        boxes_in_i = bi;
      }
    }
  }
  if(boxes_in_i<1){
    if(my_rank==0){printf("failed to find an acceptable problem size\n");}
    #ifdef USE_MPI
    MPI_Finalize();
    #endif
    exit(0);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // create the fine level...
  #ifdef USE_PERIODIC_BC
  int bc = BC_PERIODIC;
  #else
  int bc = BC_DIRICHLET;
  #endif
  level_type fine_grid;
  int ghosts=stencil_get_radius();
  create_level(&fine_grid,boxes_in_i,box_dim,ghosts,VECTORS_RESERVED,bc,my_rank,num_tasks);
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #ifdef USE_HELMHOLTZ
  double a=1.0;double b=1.0; // Helmholtz
  if(my_rank==0)fprintf(stdout,"  Creating Helmholtz (a=%f, b=%f) test problem\n",a,b);
  #else
  double a=0.0;double b=1.0; // Poisson
  if(my_rank==0)fprintf(stdout,"  Creating Poisson (a=%f, b=%f) test problem\n",a,b);
  #endif
  double h0=1.0/( (double)boxes_in_i*(double)box_dim );
  initialize_problem(&fine_grid,h0,a,b); // calculate VECTOR_ALPHA, VECTOR_BETA, and VECTOR_UTRUE
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ((a==0.0)||(dot(&fine_grid,VECTOR_ALPHA,VECTOR_ALPHA)==0.0)) && (fine_grid.boundary_condition.type == BC_PERIODIC) ){
    // Poisson w/ periodic BC's...
    // nominally, u shifted by any constant is still a valid solution.
    // However, by convention, we assume u sums to zero.
    double average_value_of_u = mean(&fine_grid,VECTOR_UTRUE);
    if(my_rank==0){fprintf(stdout,"  average value of u_true = %20.12e... shifting u_true to ensure it sums to zero...\n",average_value_of_u);}
    shift_vector(&fine_grid,VECTOR_UTRUE,VECTOR_UTRUE,-average_value_of_u);
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  //apply_op(&fine_grid,VECTOR_F,VECTOR_UTRUE,a,b); // by construction, f = A(u_true)
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(fine_grid.boundary_condition.type == BC_PERIODIC){
    double average_value_of_f = mean(&fine_grid,VECTOR_F);
    if(average_value_of_f!=0.0){
      if(my_rank==0){fprintf(stderr,"  WARNING... Periodic boundary conditions, but f does not sum to zero... mean(f)=%e\n",average_value_of_f);}
      //shift_vector(&fine_grid,VECTOR_F,VECTOR_F,-average_value_of_f);
    }
  }
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mg_type all_grids;
  int minCoarseDim = 1;
  rebuild_operator(&fine_grid,NULL,a,b); // i.e. calculate Dinv and lambda_max
  MGBuild(&all_grids,&fine_grid,a,b,minCoarseDim); // build the Multigrid Hierarchy
  double dtol=  0.0;double rtol=1e-10; // converged if ||b-Ax|| / ||b|| < rtol
//double dtol=1e-15;double rtol=  0.0; // converged if ||D^{-1}(b-Ax)|| < dtol
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     int     doTiming;
     int    minSolves = 10; // do at least minSolves MGSolves
  double timePerSolve = 0;
  for(doTiming=0;doTiming<=1;doTiming++){ // first pass warms up, second pass times

    #ifdef USE_HPM // IBM performance counters for BGQ...
    if(doTiming)HPM_Start("FMGSolve()");
    #endif

    #ifdef USE_MPI
    double minTime   = 30.0; // minimum time in seconds that the benchmark should run
    double startTime = MPI_Wtime();
    if(doTiming==1){
      if((minTime/timePerSolve)>minSolves)minSolves=(minTime/timePerSolve); // if one needs to do more than minSolves to run for minTime, change minSolves
    }
    #endif

    if(my_rank==0){
      if(doTiming==0){fprintf(stdout,"\n\n===== warming up by running %d solves ===============================\n",minSolves);}
                 else{fprintf(stdout,"\n\n===== running %d solves =============================================\n",minSolves);}
      fflush(stdout);
    }

    int numSolves =  0; // solves completed
    MGResetTimers(&all_grids);
    while( (numSolves<minSolves) ){
      zero_vector(all_grids.levels[0],VECTOR_U);
      #ifdef USE_FCYCLES
      FMGSolve(&all_grids,0,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
      #else
       MGSolve(&all_grids,0,VECTOR_U,VECTOR_F,a,b,dtol,rtol);
      #endif
      numSolves++;
    }

    #ifdef USE_MPI
    if(doTiming==0){
      double endTime = MPI_Wtime();
      timePerSolve = (endTime-startTime)/numSolves;
      MPI_Bcast(&timePerSolve,1,MPI_DOUBLE,0,MPI_COMM_WORLD); // after warmup, process 0 broadcasts the average time per solve (consensus)
    }
    #endif

    #ifdef USE_HPM // IBM performance counters for BGQ...
    if(doTiming)HPM_Stop("FMGSolve()");
    #endif
  }
  MGPrintTiming(&all_grids); // don't include the error check in the timing results
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(my_rank==0){fprintf(stdout,"calculating error...  ");}
  double fine_error = error(&fine_grid,VECTOR_U,VECTOR_UTRUE);
  if(my_rank==0){fprintf(stdout,"h = %22.15e  ||error|| = %22.15e\n\n",h0,fine_error);fflush(stdout);}
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MGDestroy()
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  #ifdef USE_MPI
  #ifdef USE_HPM // IBM performance counters for BGQ...
  HPM_Print();
  #endif
  MPI_Finalize();
  #endif
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  return;
}
