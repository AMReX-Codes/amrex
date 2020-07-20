module my_amr_module

  use iso_c_binding
  use amrex_amr_module
  use amrex_fort_module, only : rt => amrex_real

  use amr_data_module

  implicit none

  ! runtime parameters
  integer :: verbose    = 0

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
  character(len=:), allocatable, save :: check_file
  character(len=:), allocatable, save :: plot_file
  character(len=:), allocatable, save :: restart

  integer, allocatable, save :: stepno(:)
  integer, allocatable, save :: nsubsteps(:)

  real(rt), allocatable, save :: dt(:)

  integer, private, parameter :: ncomp = 1, nghost = 0
  
contains

  subroutine my_amr_init ()
    use bc_module, only : lo_bc, hi_bc
    type(amrex_parmparse) :: pp
    integer :: ilev

    if (.not.amrex_amrcore_initialized()) call amrex_amrcore_init()
    
    call amrex_init_virtual_functions (my_make_new_level_from_scratch, &
         &                             my_make_new_level_from_coarse,  &
         &                             my_remake_level,                &
         &                             my_clear_level,                 &
         &                             my_error_estimate)

    ! some default parameters
    allocate(character(len=3)::check_file)
    check_file = "chk"
    allocate(character(len=3)::plot_file)
    plot_file = "plt"
    allocate(character(len=0)::restart)

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
    call pp%query("v", verbose)
    call pp%query("verbose", verbose)
    call pp%query("cfl", cfl)
    call pp%query("do_reflux", do_reflux)
    call amrex_parmparse_destroy(pp)

    if (.not. amrex_is_all_periodic()) then
       lo_bc = amrex_bc_foextrap
       hi_bc = amrex_bc_foextrap
    end if

    allocate(stepno(0:amrex_max_level))
    stepno = 0

    allocate(nsubsteps(0:amrex_max_level))
    nsubsteps(0) = 1
    do ilev = 1, amrex_max_level
       nsubsteps(ilev) = amrex_ref_ratio(ilev-1)
    end do

    allocate(dt(0:amrex_max_level))
    dt = huge(1._rt)

    call amr_data_init()
  end subroutine my_amr_init


  subroutine my_amr_finalize ()
    call amr_data_finalize()
  end subroutine my_amr_finalize

  ! Make a new level from scratch and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_scratch (lev, time, pba, pdm) bind(c)
    use prob_module, only : init_prob_data
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phi(:,:,:,:)

    ba = pba
    dm = pdm

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call my_clear_level(lev)
  
    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)

   if (lev > 0 .and. do_reflux) then
      call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
   end if

    call amrex_mfiter_build(mfi, phi_new(lev))

    do while (mfi%next())
       bx = mfi%tilebox()
       phi => phi_new(lev)%dataptr(mfi)
       call init_prob_data(lev, t_new(lev), bx%lo, bx%hi, phi, lbound(phi), ubound(phi), &
            amrex_geom(lev)%dx, amrex_problo)
    end do

    call amrex_mfiter_destroy(mfi)
  end subroutine my_make_new_level_from_scratch

  ! Make a new level from coarse level and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_make_new_level_from_coarse (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillcoarsepatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm

    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm

    ba = pba
    dm = pdm

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    call fillcoarsepatch(lev, time, phi_new(lev))
  end subroutine my_make_new_level_from_coarse

  ! Remake a level from current and coarse elvels and put the data in phi_new.
  ! Note tha phi_old contains no valid data after this.
  subroutine my_remake_level (lev, time, pba, pdm) bind(c)
    use fillpatch_module, only : fillpatch
    integer, intent(in), value :: lev
    real(amrex_real), intent(in), value :: time
    type(c_ptr), intent(in), value :: pba, pdm
    
    type(amrex_boxarray) :: ba
    type(amrex_distromap) :: dm
    type(amrex_multifab) :: new_phi_new

    ba = pba
    dm = pdm

    call amrex_multifab_build(new_phi_new, ba, dm, ncomp, 0)
    call fillpatch(lev, time, new_phi_new)

    call my_clear_level(lev)

    t_new(lev) = time
    t_old(lev) = time - 1.e200_amrex_real

    call amrex_multifab_build(phi_new(lev), ba, dm, ncomp, nghost)
    call amrex_multifab_build(phi_old(lev), ba, dm, ncomp, nghost)
    if (lev > 0 .and. do_reflux) then
       call amrex_fluxregister_build(flux_reg(lev), ba, dm, amrex_ref_ratio(lev-1), lev, ncomp)
    end if

    call phi_new(lev)%copy(new_phi_new, 1, 1, ncomp, 0)

    call amrex_multifab_destroy(new_phi_new)
  end subroutine my_remake_level

  subroutine my_clear_level (lev) bind(c)
    integer, intent(in), value :: lev
    call amrex_multifab_destroy(phi_new(lev))
    call amrex_multifab_destroy(phi_old(lev))
    call amrex_fluxregister_destroy(flux_reg(lev))
  end subroutine my_clear_level

  subroutine my_error_estimate (lev, cp, t, settag, cleartag) bind(c)
    use tagging_module, only : tag_phi_error
    integer, intent(in), value :: lev
    type(c_ptr), intent(in), value :: cp
    real(amrex_real), intent(in), value :: t
    character(kind=c_char), intent(in), value :: settag, cleartag

    real(amrex_real), allocatable, save :: phierr(:)
    type(amrex_parmparse) :: pp
    type(amrex_tagboxarray) :: tag
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx
    real(amrex_real), contiguous, pointer :: phiarr(:,:,:,:)
    character(kind=c_char), contiguous, pointer :: tagarr(:,:,:,:)

    if (.not.allocated(phierr)) then
       call amrex_parmparse_build(pp, "myamr")
       call pp%getarr("phierr", phierr)
       call amrex_parmparse_destroy(pp)
    end if

    tag = cp

    !$omp parallel private(mfi, bx, phiarr, tagarr)
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()
       phiarr => phi_new(lev)%dataptr(mfi)
       tagarr => tag%dataptr(mfi)
       call tag_phi_error(lev, t, bx%lo, bx%hi, &
            phiarr, lbound(phiarr), ubound(phiarr), &
            tagarr, lbound(tagarr), ubound(tagarr), &
            phierr(lev+1), settag, cleartag)  ! +1 because level starts with 0, but phierr starts with 1
    end do
    call amrex_mfiter_destroy(mfi)
    !$omp end parallel

  end subroutine my_error_estimate
  
end module my_amr_module
