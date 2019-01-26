//------------------------------------------------------------------------------------------------------------------------------
// Samuel Williams
// SWWilliams@lbl.gov
// Lawrence Berkeley National Lab
//------------------------------------------------------------------------------------------------------------------------------
#ifdef _OPENMP
// getTime in OpenMP is now defined as a preprocessor macro
//#include "./timers/omp.c"
#elif USE_MPI
// getTime in MPI is now defined as a preprocessor macro
//#include "./timers/mpi.c"
#else
#include "./timers/x86.c"
#endif
