module my_amr_module

  use amrex_module
  use amrex_famrcore_module

  implicit none

  private

  ! runtime parameters
  integer :: max_step   = huge(1)
  integer :: max_level  = 0
  integer :: regrid_int = 2
  integer :: check_int  = -1
  integer :: plot_int   = -1
  !
  logical :: do_reflux  = .true.
  !
  real(amrex_real) :: stop_time  = huge(1._amrex_real)
  real(amrex_real) :: cfl        = 0.7_amrex_real
  !
  character(len=127) :: check_file = "chk"
  character(len=127) :: plot_file  = "plt"
  character(len=127) :: restart    = ""  

  integer, allocatable :: istep(:)
  integer, allocatable :: nsubsteps(:)
  
  public :: my_amr_init, my_amr_finalize

contains

  subroutine my_amr_init ()
    type(amrex_parmparse) :: pp

    if (.not.amrex_famrcore_initialized()) call amrex_famrcore_init()

    ! Read parameters
    call amrex_parmparse_build(pp)
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters amr.*
    call amrex_parmparse_build(pp, "amr")
    call pp%query("max_level", max_level)
    call pp%query("regrid_int", regrid_int)
    call pp%query("check_int", check_int)
    call pp%query("plot_int", plot_int)
    call pp%query("check_file", check_file)
    call pp%query("plot_file", plot_file)
    call pp%query("restart", restart)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters myamr.*
    call amrex_parmparse_build(pp, "myamr")
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call amrex_parmparse_destroy(pp)

  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    
  end subroutine my_amr_finalize

end module my_amr_module
