
module fab_module

  use bl_space_module, only : ndims => bl_num_dims

  implicit none

  private

  type, public :: Fab
     integer, dimension(3) :: lo = 1
     integer, dimension(3) :: hi = 1
     integer               :: nc = 1
     double precision, pointer, dimension(:,:,:,:) :: p => null()
  end type Fab

end module fab_module
