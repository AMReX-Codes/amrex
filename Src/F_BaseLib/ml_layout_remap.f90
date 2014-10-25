module ml_layout_remap_module

  use layout_module
  use multifab_module

  implicit none

  private

  public :: ml_layout_remap

contains

  subroutine ml_layout_remap(la, phi, nlevs, rr)
    type(layout), intent(inout) :: la(:)
    type(multifab), intent(inout) :: phi(:)
    integer, intent(in) :: nlevs, rr(:,:)

    
    
  end subroutine ml_layout_remap

end module ml_layout_remap_module
