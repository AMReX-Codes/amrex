module fillpatch_module

  use iso_c_binding
  use amrex_amr_module

  implicit none

  private

  public :: fillpatch

contains

  ! Fill phi with data from phi_old and phi_new of current level and one level below.
  subroutine fillpatch (lev, time, phi)
    use my_amr_module, only : t_old, t_new, phi_old, phi_new
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1  ! for this test code
    real(amrex_real) :: teps
    type(amrex_multifab), allocatable :: c_mf(:), f_mf(:)
    real(amrex_real), allocatable :: c_time(:), f_time(:)
    type(amrex_physbc) :: c_pbc, f_pbc
    integer :: lo_bc(amrex_spacedim,src_comp+num_comp-1)
    integer :: hi_bc(amrex_spacedim,src_comp+num_comp-1)

    teps = 1.e-3_amrex_real * (t_new(lev) - t_old(lev))
    if (abs(time-t_new(lev)) .lt. teps) then
       allocate(f_mf(1))
       allocate(f_time(1))
       f_mf(1) = phi_new(lev)
       f_time(1) = t_new(lev)
    else if (abs(time-t_old(lev)) .lt. teps) then
       allocate(f_mf(1))
       allocate(f_time(1))
       f_mf(1) = phi_old(lev)
       f_time(1) = t_old(lev)
    else
       allocate(f_mf(2))
       allocate(f_time(2))
       f_mf(1) = phi_old(lev)
       f_mf(2) = phi_new(lev)
       f_time(1) = t_old(lev)
       f_time(2) = t_new(lev)
    end if
 
    if (lev .eq. 0) then
       call amrex_fillpatch(phi, t_old(lev), phi_old(lev), &
            &                    t_new(lev), phi_new(lev), &
            &               amrex_geom(lev), fill_physbc , &
            time, src_comp, dst_domp, num_comp)
    else
       call amrex_fillpatch(phi, t_old(lev-1), phi_old(lev-1), &
            &                    t_new(lev-1), phi_new(lev-1), &
            &               amrex_geom(lev-1), fill_physbc   , &
            &                    t_old(lev  ), phi_old(lev  ), &
            &                    t_new(lev  ), phi_new(lev  ), &
            &               amrex_geom(lev  ), fill_physbc   , &
            time, src_comp, dst_domp, num_comp)


       teps = 1.e-3_amrex_real * (t_new(lev-1) - t_old(lev-1))
       if (abs(time-t_new(lev-1)) .lt. teps) then
          allocate(c_mf(1))
          allocate(c_time(1))
          c_mf(1) = phi_new(lev-1)
          c_time(1) = t_new(lev-1)
       else if (abs(time-t_old(lev-1)) .lt. teps) then
          allocate(c_mf(1))
          allocate(c_time(1))
          c_mf(1) = phi_old(lev-1)
          c_time(1) = t_old(lev-1)
       else
          allocate(c_mf(2))
          allocate(c_time(2))
          c_mf(1) = phi_old(lev-1)
          c_mf(2) = phi_new(lev-1)
          c_time(1) = t_old(lev-1)
          c_time(2) = t_new(lev-1)
       end if

       call amrex_physbc_build(c_pbc, fill_physbc)
       call amrex_physbc_build(f_pbc, fill_physbc)

       ! periodic bc
       lo_bc = amrex_bc_int_dir
       hi_bc = amrex_bc_int_dir

       call amrex_fillpatch(phi, time, &
            c_mf, c_time, f_mf, f_time, &
            src_comp, dst_comp, num_comp, &
            amrex_geom(lev-1), amrex_geom(lev), c_pbc, f_pbc, &
            amrex_ref_ratio(lev-1), amrex_interp_cell_cons, lo_bc, hi_bc)

       call amrex_physbc_destroy(f_pbc)
       call amrex_physbc_destroy(c_pbc)
    end if
  end subroutine fillpatch


  subroutine fill_physbc (pmf, scomp, ncomp, time) bind(c)
    type(c_ptr), value :: pmf
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time
    ! In this test problem, we only have periodic boundaries.
    ! So there is noting to do.
    
!    type(amrex_multifab) :: mf
!    mf = pmf
    
  end subroutine fill_physbc

end module fillpatch_module
