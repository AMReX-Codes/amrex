#include "HYPRE_config.h"
module amrex_hypre_fort_module
  use iso_c_binding, only : c_long_long, c_int
  implicit none
#ifdef HYPRE_BIGINT
  integer, parameter :: hypre_int = c_long_long
#else
  integer, parameter :: hypre_int = c_int
#endif
end module amrex_hypre_fort_module
