//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#include <stdint.h>
#include <mpi.h>
double getTime(){
  return(MPI_Wtime()); // timers are in units of seconds; no conversion is necessary
}
