#include <AMReX_EBStruct.H>

module amrex_ebstruct_module

  use amrex_fort_module, only : dim=>bl_spacedim
  implicit none

  type, bind(c) :: cutface
     integer :: cellHi, cellLo
     integer :: parent
  end type cutface

  type, bind(c) :: cutcell
     integer :: Nnbr(0:1,0:dim-1)
     integer :: nbr(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: faceID(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: parent
  end type cutcell

end module amrex_ebstruct_module
