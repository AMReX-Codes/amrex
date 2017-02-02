
module amrex_box_module

  use amrex_fort_module, only : ndims => amrex_spacedim;

  implicit none

  private

  type, public :: amrex_box
     integer, dimension(3) :: lo = 1
     integer, dimension(3) :: hi = 1
  end type amrex_box

  interface amrex_box
     module procedure amrex_box_build
  end interface amrex_box

contains

  function amrex_box_build(lo, hi) result(bx)
    integer, intent(in) :: lo(*), hi(*)
    type(amrex_box) :: bx
    bx%lo(1:ndims) = lo(1:ndims)
    bx%hi(1:ndims) = hi(1:ndims)
  end function amrex_box_build

end module amrex_box_module

