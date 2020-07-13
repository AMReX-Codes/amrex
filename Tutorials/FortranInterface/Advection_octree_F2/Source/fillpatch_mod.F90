module fillpatch_module

  use iso_c_binding
  use amrex_amr_module

  implicit none

  private

  public :: fillpatch, fillcoarsepatch

contains

  ! Fill phi with data from phi_old and phi_new of current level and one level below.
  subroutine fillpatch (lev, time, phi)
    use amr_data_module, only : t_old, t_new, phi_old, phi_new
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi
    
    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1  ! for this test code

    if (lev .eq. 0) then
       call amrex_fillpatch(phi, t_old(lev), phi_old(lev), &
            &                    t_new(lev), phi_new(lev), &
            &               amrex_geom(lev), fill_physbc , &
            &               time, src_comp, dst_comp, num_comp)
    else
       call amrex_fillpatch(phi, t_old(lev-1), phi_old(lev-1), &
            &                    t_new(lev-1), phi_new(lev-1), &
            &               amrex_geom(lev-1), fill_physbc   , &
            &                    t_old(lev  ), phi_old(lev  ), &
            &                    t_new(lev  ), phi_new(lev  ), &
            &               amrex_geom(lev  ), fill_physbc   , &
            &               time, src_comp, dst_comp, num_comp, &
            &               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
            &               lo_bc, hi_bc, pre_interp, post_interp)
       ! see amrex_interpolater_module for a list of interpolaters
    end if
  end subroutine fillpatch

  subroutine fillcoarsepatch (lev, time, phi)
    use amr_data_module, only : t_old, t_new, phi_old, phi_new
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1  ! for this test code
    
    call amrex_fillcoarsepatch(phi, t_old(lev-1), phi_old(lev-1),  &
         &                          t_new(lev-1), phi_new(lev-1),  &
         &                     amrex_geom(lev-1),    fill_physbc,  &
         &                     amrex_geom(lev  ),    fill_physbc,  &
         &                     time, src_comp, dst_comp, num_comp, &
         &                     amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
         &                     lo_bc, hi_bc, pre_interp, post_interp)
       ! see amrex_interpolater_module for a list of interpolaters
  end subroutine fillcoarsepatch

  subroutine fill_physbc (pmf, scomp, ncomp, time, pgeom) bind(c)
    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time
    ! In this test problem, we only have periodic boundaries.
    ! So there is noting to do.
    
!    type(amrex_multifab) :: mf
!    mf = pmf
    
  end subroutine fill_physbc

  subroutine pre_interp (lo, hi, d, dlo, dhi, nd, icomp, ncomp) bind(c)
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    integer(c_int), intent(in), value :: nd, icomp, ncomp
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nd)

    ! one might modify d(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),icomp:icomp+ncomp-1)
    !
    ! In 2d, lo(3) = hi(3) = dlo(3) = dhi(3) = 0
  end subroutine pre_interp

  subroutine post_interp (lo, hi, d, dlo, dhi, nd, icomp, ncomp) bind(c)
    integer(c_int), intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    integer(c_int), intent(in), value :: nd, icomp, ncomp
    real(amrex_real), intent(inout) :: d(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),nd)

    ! one might modify d(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),icomp:icomp+ncomp-1)
    !
    ! In 2d, lo(3) = hi(3) = dlo(3) = dhi(3) = 0
  end subroutine post_interp

end module fillpatch_module
