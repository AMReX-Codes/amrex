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

    integer, parameter :: src_comp=0, dst_comp=0, num_comp=1  ! for this test code
    real(amrex_real) :: teps
    type(amrex_multifab), allocatable :: crse_mf(:), curr_mf(:)
    real(amrex_real), allocatable :: crse_time(:), curr_time(:)
    type(amrex_physbc) :: crse_pbc, curr_pbc

    if (lev .eq. 0) then
       teps = 1.e-3_amrex_real * (t_new(lev) - t_old(lev))
       if (abs(time-t_new(lev)) .lt. teps) then
          allocate(curr_mf(1))
          allocate(curr_time(1))
          curr_mf(1) = phi_new(lev)
          curr_time(1) = t_new(lev)
       else if (abs(time-t_old(lev)) .lt. teps) then
          allocate(curr_mf(1))
          allocate(curr_time(1))
          curr_mf(1) = phi_old(lev)
          curr_time(1) = t_old(lev)
       else
          allocate(curr_mf(2))
          allocate(curr_time(2))
          curr_mf(1) = phi_old(lev)
          curr_mf(2) = phi_new(lev)
          curr_time(1) = t_old(lev)
          curr_time(2) = t_new(lev)
       end if
 
       call amrex_physbc_build(curr_pbc, fill_physbc)
      
       call amrex_fillpatch(phi, time, curr_mf, curr_time, src_comp, dst_comp, num_comp, &
            &               amrex_geom(lev), curr_pbc)


       call amrex_physbc_destroy(curr_pbc)
    else
!       xxxxxx
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
