#include <petsc/finclude/petscsys.h>

module amrex_petsc_fort_module
  implicit none
  PetscInt, parameter, private :: a_petsc_int=0
  integer, parameter, public :: petsc_int = kind(a_petsc_int)
end module amrex_petsc_fort_module
