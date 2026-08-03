module fillpatch_module

  use iso_c_binding
  use amrex_amr_module

  implicit none

  private

  public :: fillpatch, fillcoarsepatch, fillpatch_stage
  public :: store_rk3_fillpatcher, store_rk4_fillpatcher
  public :: reset_fillpatcher, fill_physbc

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
            &               lo_bc, hi_bc)
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
         &                     lo_bc, hi_bc)
       ! see amrex_interpolater_module for a list of interpolaters
  end subroutine fillcoarsepatch

  subroutine fillpatch_stage (lev, time, phi, rk_order, stage, iteration, ncycle)
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, fillpatcher
    use bc_module, only : lo_bc, hi_bc
    integer, intent(in) :: lev, rk_order, stage, iteration, ncycle
    real(amrex_real), intent(in) :: time
    type(amrex_multifab), intent(inout) :: phi

    integer, parameter :: src_comp=1, dst_comp=1, num_comp=1

    if (lev .eq. 0) then
       call phi%fill_boundary(amrex_geom(lev))
       call fill_physbc(phi%p, src_comp, num_comp, time, amrex_geom(lev)%p)
    else if (rk_order .ge. 3) then
       call fillpatcher(lev)%fill_rk(stage, iteration, ncycle, phi, time, &
            amrex_geom(lev-1), fill_physbc, amrex_geom(lev), fill_physbc, &
            lo_bc, hi_bc)
    else
       call amrex_fillpatch(phi, t_old(lev-1), phi_old(lev-1), &
            &                    t_new(lev-1), phi_new(lev-1), &
            &               amrex_geom(lev-1), fill_physbc,    &
            &                    time, phi,                    &
            &                    time, phi,                    &
            &               amrex_geom(lev), fill_physbc,      &
            &               time, src_comp, dst_comp, num_comp,&
            &               amrex_ref_ratio(lev-1), amrex_interp_cell_cons, &
            &               lo_bc, hi_bc)
    end if
  end subroutine fillpatch_stage

  ! Store coarse AMR level RK stage data into the FillPatcher on the next
  ! finer level so that the fine level can fill its ghost cells at the
  ! coarse/fine boundary during RK substeps (via fillpatch_stage -> fill_rk).
  !
  ! build_fillpatcher (not ensure_fillpatcher) is called here because this is
  ! the only place the FillPatcher for lev+1 is created. It is not built
  ! anywhere during fine-level advance; we must build it now. The FillPatcher
  ! is destroyed later by reset_fillpatcher after the fine level catches up.
  subroutine store_rk3_fillpatcher (lev, time, dt, s_old, rk, nghost)
    use amr_data_module, only : fillpatcher
    integer, intent(in) :: lev, nghost
    real(amrex_real), intent(in) :: time, dt
    type(amrex_multifab), intent(in) :: s_old, rk(3)

    if (lev < amrex_get_finest_level()) then
       call build_fillpatcher(lev+1, nghost)
       call fillpatcher(lev+1)%store_rk3(time, dt, s_old, rk)
    end if
  end subroutine store_rk3_fillpatcher

  subroutine store_rk4_fillpatcher (lev, time, dt, s_old, rk, nghost)
    use amr_data_module, only : fillpatcher
    integer, intent(in) :: lev, nghost
    real(amrex_real), intent(in) :: time, dt
    type(amrex_multifab), intent(in) :: s_old, rk(4)

    if (lev < amrex_get_finest_level()) then
       call build_fillpatcher(lev+1, nghost)
       call fillpatcher(lev+1)%store_rk4(time, dt, s_old, rk)
    end if
  end subroutine store_rk4_fillpatcher

  subroutine reset_fillpatcher (lev)
    use amr_data_module, only : fillpatcher
    integer, intent(in) :: lev
    if (lev >= 0 .and. lev <= amrex_max_level) then
       call amrex_fillpatcher_destroy(fillpatcher(lev))
    end if
  end subroutine reset_fillpatcher

  subroutine ensure_fillpatcher (lev, nghost)
    use amr_data_module, only : fillpatcher
    integer, intent(in) :: lev, nghost
    if (.not. c_associated(fillpatcher(lev)%p)) then
       call build_fillpatcher(lev, nghost)
    end if
  end subroutine ensure_fillpatcher

  subroutine build_fillpatcher (lev, nghost)
    use amr_data_module, only : phi_new, fillpatcher
    integer, intent(in) :: lev, nghost
    integer, parameter :: num_comp=1

    call amrex_fillpatcher_build(fillpatcher(lev), &
         phi_new(lev)%ba,   phi_new(lev)%dm,   amrex_geom(lev), &
         phi_new(lev-1)%ba, phi_new(lev-1)%dm, amrex_geom(lev-1), &
         nghost, num_comp, amrex_interp_cell_cons)
  end subroutine build_fillpatcher

  subroutine fill_physbc (pmf, scomp, ncomp, time, pgeom) bind(c)
    use amrex_geometry_module, only : amrex_is_all_periodic
    use amrex_filcc_module, only : amrex_filcc
    use bc_module, only : lo_bc, hi_bc
    type(c_ptr), value :: pmf, pgeom
    integer(c_int), value :: scomp, ncomp
    real(amrex_real), value :: time

    type(amrex_geometry) :: geom
    type(amrex_multifab) :: mf
    type(amrex_mfiter) :: mfi
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: p
    integer :: plo(4), phi(4)

    if (.not. amrex_is_all_periodic()) then
       geom = pgeom
       mf = pmf

       !$omp parallel private(mfi,p,plo,phi)
       call amrex_mfiter_build(mfi, mf, tiling=.false.)
       do while(mfi%next())
          p => mf%dataptr(mfi)
          if (.not. geom%domain%contains(p)) then ! part of this box is outside the domain
             plo = lbound(p)
             phi = ubound(p)
             call amrex_filcc(p, plo, phi,         & ! fortran array and bounds
                  geom%domain%lo, geom%domain%hi,  & ! index extent of whole problem domain
                  geom%dx,                         & ! cell size in real
                  geom%get_physical_location(plo), & ! physical location of lower left corner
                  lo_bc, hi_bc)                      ! bc types for each component

             ! amrex_filcc doesn't fill EXT_DIR/EXT_DIR_CC (see amrex_bc_types_module for a list of bc types
             ! In that case, the user needs to fill it.
          end if
       end do
       !$omp end parallel

    end if

  end subroutine fill_physbc

end module fillpatch_module
