module initdata_module

  use my_amr_module

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then
       call amrex_init_from_scratch(0.0_amrex_real)
!xxx       call average_down()
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

end module initdata_module
