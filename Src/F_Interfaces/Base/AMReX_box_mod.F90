
module box_module

  use bl_space_module, only : ndims => bl_num_dims

  implicit none

  private

  type, public :: Box
     integer, dimension(3) :: lo = 1
     integer, dimension(3) :: hi = 1
  end type Box

  interface Box
     module procedure build_box
  end interface Box

contains

  function build_box(lo, hi) result(bx)
    integer, intent(in) :: lo(:), hi(:)
    type(Box) :: bx
    bx%lo(1:ndims) = lo(1:ndims)
    bx%hi(1:ndims) = hi(1:ndims)
  end function build_box

end module box_module

