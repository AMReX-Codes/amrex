
module amrex_box_module

  use amrex_fort_module, only : ndims => amrex_spacedim;

  implicit none

  private

  type, public :: amrex_box
     integer, dimension(3) :: lo = 1
     integer, dimension(3) :: hi = 1
  end type amrex_box

  public :: amrex_print, amrex_allprint

  interface amrex_box
     module procedure amrex_box_build
  end interface amrex_box

  interface amrex_print
     module procedure amrex_box_print
  end interface amrex_print

  interface amrex_allprint
     module procedure amrex_box_allprint
  end interface amrex_allprint

  interface
     subroutine amrex_fi_print_box(lo,hi,all) bind(c)
       use iso_c_binding
       integer(c_int), intent(in) :: lo(*), hi(*)
       integer(c_int), value :: all
     end subroutine amrex_fi_print_box
  end interface

contains

  function amrex_box_build (lo, hi) result(bx)
    integer, intent(in) :: lo(*), hi(*)
    type(amrex_box) :: bx
    bx%lo(1:ndims) = lo(1:ndims)
    bx%hi(1:ndims) = hi(1:ndims)
  end function amrex_box_build

  subroutine amrex_box_print (bx)
    type(amrex_box), intent(in) :: bx
    call amrex_fi_print_box(bx%lo, bx%hi, 0)
  end subroutine amrex_box_print

  subroutine amrex_box_allprint (bx)
    type(amrex_box), intent(in) :: bx
    call amrex_fi_print_box(bx%lo, bx%hi, 1)
  end subroutine amrex_box_allprint

end module amrex_box_module

