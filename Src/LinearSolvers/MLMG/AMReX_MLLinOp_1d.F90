
module mllinop_1d_module
  use amrex_fort_module, only : amrex_real
  implicit none

  private
  public :: amrex_mllinop_grad

contains
  
  subroutine amrex_mllinop_grad () bind(c,name='amrex_mllinop_grad');

  end subroutine amrex_mllinop_grad

end module mllinop_1d_module
