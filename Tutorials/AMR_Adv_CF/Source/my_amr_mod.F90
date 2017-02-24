module my_amr_module

  use amrex_module
  use amrex_famrcore_module

  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! runtime parameters
  integer :: max_step   = huge(1)
  integer :: regrid_int = 2
  integer :: check_int  = -1
  integer :: plot_int   = -1
  !
  logical :: do_reflux  = .true.
  !
  real(rt) :: stop_time  = huge(1._rt)
  real(rt) :: cfl        = 0.7_rt
  !
  character(len=127) :: check_file = "chk"
  character(len=127) :: plot_file  = "plt"
  character(len=127) :: restart    = ""  

  integer, allocatable :: ref_ratio(:)
  integer, allocatable :: istep(:)
  integer, allocatable :: nsubsteps(:)
  integer, allocatable :: last_regrid_step(:)

  real(rt), allocatable :: t_new(:)
  real(rt), allocatable :: t_old(:)
  real(rt), allocatable :: dt(:)

  type(amrex_multifab), allocatable :: phi_new(:)
  type(amrex_multifab), allocatable :: phi_old(:)

  ! type(amrex_fluxregister), allocatable :: flux_reg(:)
  
contains

  subroutine my_amr_init ()
    type(amrex_parmparse) :: pp
    integer :: ilev

    if (.not.amrex_famrcore_initialized()) call amrex_famrcore_init()
    
    call amrex_init_virtual_functions (c_funloc(my_make_new_level_from_scratch))

    ! Read parameters
    call amrex_parmparse_build(pp)
    call pp%query("max_step", max_step)
    call pp%query("stop_time", stop_time)
    call amrex_parmparse_destroy(pp)
    
    ! Parameters amr.*
    call amrex_parmparse_build(pp, "amr")
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

    allocate(istep(0:amrex_max_level))
    istep = 0

    allocate(nsubsteps(0:amrex_max_level))
    nsubsteps(0) = 1
    do ilev = 1, amrex_max_level
       nsubsteps(ilev) = amrex_ref_ratio(ilev-1)
    end do

    allocate(last_regrid_step(0:amrex_max_level))
    last_regrid_step = 0

    allocate(t_new(0:amrex_max_level))
    t_new = 0.0_rt

    allocate(t_old(0:amrex_max_level))
    t_old = -1.0e100_rt

    allocate(dt(0:amrex_max_level))
    dt = 1.e100

    allocate(phi_new(0:amrex_max_level))
    allocate(phi_old(0:amrex_max_level))

    ! allocate(flux_reg(0:amrex_max_level))

  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    integer :: lev
    do lev = 0, amrex_max_level
       call amrex_multifab_destroy(phi_new(lev))
       call amrex_multifab_destroy(phi_old(lev))
       ! call amrex_fluxregister_destroy(flux_reg(lev))
    end do
  end subroutine my_amr_finalize


  subroutine make_new_level (lev, time, ba, dm)
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    
    integer, parameter :: ncomp = 1, nghost = 0

    call amrex_install_level(lev, ba, dm)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real
  
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)

!    if (lev > 0 .and. do_reflux) then
!       call amrex_fluxregister_destroy(flux_reg(lev))
!       call amrex_fluxregister_build(flux_reg(lev), ba, dm, ref_ratio(lev-1), lev, ncomp)
!    end if
  end subroutine make_new_level

  subroutine my_make_new_level_from_scratch (lev, time, ba, dm) bind(c)
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), value :: ba, dm
    print *, 'in my_make_new_level_from_scratch'
  end subroutine my_make_new_level_from_scratch
  
end module my_amr_module
