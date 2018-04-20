#include <stdio.h>
#include "AMReX_PXStuff.H"

namespace amrex
{

  const char *PXE_ErrorCodeName[17] = { "No Error",
                                        "Memory",
                                        "Bad Input",
                                        "Nonphysical",
                                        "Read/Write",
                                        "Grid",
                                        "Search Not Found",
                                        "No Update",
                                        "Parallel",
                                        "Code Flow",
                                        "System",
                                        "Dynamic Library",
                                        "Not Converged",
                                        "Viz",
                                        "Lapack",
                                        "Hard Exit",
                                        "CGNS"};

  void PXErrorReport( const char *file, const int line, const char *call, const int ierr){

/*   if (ierr == PX_NO_ERROR) */
/*     return PX_NO_ERROR; */

    printf("Error %d (%s) has occured.\n File : %s  Line : %d\n Call : %s\n", ierr,PXE_ErrorCodeName[-ierr\
             ], file, line, call);
    //  printf("Error %d has occured.\n   File : %s   Line : %d\n   Call : %s\n", ierr, file, line, call);
    fflush(stdout);
  }

}
