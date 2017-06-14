#include <Node.H>

module amrex_ebstruct_module

  use amrex_fort_module, only : dim=>bl_spacedim
  implicit none

  type, bind(c) :: cutface
     integer :: cellHi, cellLo
     integer :: ebFaceID
  end type cutface

  type, bind(c) :: cutcell
     integer :: Nnbr(0:1,0:dim-1)
     integer :: nbr(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: faceID(0:NCELLMAX-1,0:1,0:dim-1)
     integer :: ebCellID
  end type cutcell

  type, bind(c) :: fnode
     integer :: nFaces
     integer :: iv(0:dim-1)
     type(cutface) :: faces(0:NFACEMAX-1)
  end type fnode

  type, bind(c) :: cnode
     integer :: nCells
     integer :: iv(0:dim-1)
     type(cutcell) :: cells(0:NCELLMAX-1)
  end type cnode

end module amrex_ebstruct_module
