//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifndef TIMER_H
#define TIMER_H

  #include<stdint.h>

  #ifdef _OPENMP
    #include <omp.h>
    #define getTime() (omp_get_wtime())

  #elif USE_MPI
    #include <mpi.h>
    #define getTime() (MPI_Wtime())

  #else
    // user must provide a function getTime and include it in timers.c
    // if calibration is necesary, then the user must #define CALIBRATE_TIMER
    double getTime();
  #endif

#endif
