module ml_remap_module

  use layout_module
  use ml_layout_module
  use multifab_module
  use ml_optimization_module

  implicit none

  private

  public :: ml_remap

contains

  subroutine ml_remap(la, phi, nlevs, rr)
    type(layout), intent(inout) :: la(:)
    type(multifab), intent(inout) :: phi(:)
    integer, intent(in) :: nlevs, rr(:,:)

    type(layout) :: la_new(nlevs)

    call optimiza_ml_layout(la_new, la, nlevs, rr)
    
    ! test if la_new(n) == la(n)
    ! if ne, we need to build phi_new, and copy, and install phi_new as phi

  end subroutine ml_remap

end module ml_remap_module
