#include <iostream>
#include <mpi.h>
#include <BoxLib/Profiler.H>

{{fnall fn_name MPI_Type_count MPI_Wtime }}
  BL_PROFILE_TIMER( bltimer, "{{fn_name}}()" );	
  BL_PROFILE_START( bltimer );
  {{callfn}}
  BL_PROFILE_STOP( bltimer );
{{endfnall}}
