module initdata_module

  use my_amr_module

  implicit none

  private

  public :: initdata

contains

  subroutine initdata ()
    if (len_trim(restart) .eq. 0) then
       call init_from_scratch();
    else
       call amrex_abort("init from checkpoint not implemented yet")
    end if
  end subroutine initdata

  subroutine init_from_scratch ()
    type(amrex_boxarray) :: ba0
    type(amrex_distromap) :: dm0
    call amrex_make_base_grids(ba0)
    call amrex_distromap_build(dm0, ba0)
    call amrex_print(dm0)
  end subroutine init_from_scratch

end module initdata_module
