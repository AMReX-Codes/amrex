module evolve_module

  use amrex_amr_module

  use amrex_particlecontainer_module, only : amrex_particle
  
  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()
    use my_amr_module, only : stepno, max_step, stop_time, dt, plot_int
    use amr_data_module, only : phi_old, phi_new, t_new
    use compute_dt_module, only : compute_dt
    use plotfile_module, only : writeplotfile
    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, substep, finest_level

    cur_time = t_new(0)
    last_plot_file_step = 0;
    
    do step = stepno(0), max_step-1
       if (cur_time .ge. stop_time) exit

       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       call compute_dt()

       lev = 0
       substep = 1
       call timestep(lev, cur_time, substep)
 
       cur_time = cur_time + dt(0)

       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dt(0)
       end if

       ! sync up time to avoid roundoff errors
       finest_level = amrex_get_finest_level()
       do lev = 0, finest_level
          t_new(lev) = cur_time
       end do

       if (plot_int .gt. 0 .and. mod(step+1,plot_int) .eq. 0) then
          last_plot_file_step = step+1
          call writeplotfile()
       end if

       if (cur_time .ge. stop_time - 1.e-6_amrex_real*dt(0)) exit
    end do

    if (plot_int .gt. 0 .and. stepno(0) .gt. last_plot_file_step) then
       call writeplotfile()
    end if

  end subroutine evolve

  recursive subroutine timestep (lev, time, substep)
    use my_amr_module, only : regrid_int, stepno, nsubsteps, dt, do_reflux
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg, pc
    use averagedown_module, only : averagedownto
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    
    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep
    integer :: redistribute_ngrow
    
    if (regrid_int .gt. 0) then
       if (.not.allocated(last_regrid_step)) then
          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
       end if

       ! regrid doesn't change the base level, so we don't regrid on amrex_max_level
       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then
          if (mod(stepno(lev), regrid_int) .eq. 0) then

             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time) ! the finest level may change during regrid
             finest_level = amrex_get_finest_level()

             do k = lev, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
             end do

             call pc%redistribute(lev)
          end if
       end if
    end if

    stepno(lev) = stepno(lev)+1

    ! We need to update t_old(lev) and t_new(lev) before advance is called because of fillpath.
    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    ! swap phi_new(lev) and phi_old(lev) so they are consistent with t_new(lev) and t_old(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))

    call advance(lev, time, dt(lev), stepno(lev), substep, nsubsteps(lev))

    if (lev .lt. amrex_get_finest_level()) then
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep)
       end do

       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real)
       end if

       call averagedownto(lev)
    end if

    if ((substep < nsubsteps(lev)) .or. lev == 0) then
       if (lev == 0) then
          redistribute_ngrow = 0
       else
          redistribute_ngrow = substep
       end if
       call pc%redistribute(lev, amrex_get_finest_level(), redistribute_ngrow)
    end if
  
  end subroutine timestep

  ! update phi_new(lev)
  subroutine advance (lev, time, dt, step, substep, nsub)
    use my_amr_module, only : verbose, do_reflux
    use amr_data_module, only : phi_new, flux_reg, pc
    use face_velocity_module, only : get_face_velocity
    use advect_module, only : advect, advect_particles
    use fillpatch_module, only : fillpatch
    integer, intent(in) :: lev, step, substep, nsub
    real(amrex_real), intent(in) :: time, dt

    integer, parameter :: ngrow = 3
    integer :: ncomp, idim
    logical :: nodal(3)
    type(amrex_multifab) :: phiborder
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout,pux,puy,puz,pfx,pfy,pfz, &
         pf, pfab
    type(amrex_fab) :: uface(amrex_spacedim)
    type(amrex_fab) ::  flux(amrex_spacedim)
    type(amrex_multifab) :: fluxes(amrex_spacedim)
    type(amrex_particle), pointer :: particles(:)
    
    if (verbose .gt. 0 .and. amrex_parallel_ioprocessor()) then
       write(*,'(A, 1X, I0, 1X, A, 1X, I0, A, 1X, G0)') &
            "[Level", lev, "step", step, "] ADVANCE with dt =", dt
    end if

    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)

    call fillpatch(lev, time, phiborder)

    !$omp parallel private(mfi,bx,tbx,pin,pout,pux,puy,puz,pfx,pfy,pfz,pf,pfab,uface,flux,particles)
    do idim = 1, amrex_spacedim
       call uface(idim)%reset_omp_private()
       call flux(idim)%reset_omp_private()
    end do
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()

       pin  => phiborder%dataptr(mfi)
       pout => phi_new(lev)%dataptr(mfi)

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
          call tbx%grow(substep)
          call uface(idim)%resize(tbx,1)
       end do

       pux => uface(1)%dataptr()
       pfx =>  flux(1)%dataptr()
       puy => uface(2)%dataptr()
       pfy =>  flux(2)%dataptr()
#if BL_SPACEDIM == 3
       puz => uface(3)%dataptr()
       pfz =>  flux(3)%dataptr()
#endif

       call get_face_velocity(time+0.5_amrex_real*dt, &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            amrex_geom(lev)%dx, amrex_problo)

       call advect(time, bx%lo, bx%hi, &
            pin, lbound(pin), ubound(pin), &
            pout,lbound(pout),ubound(pout), &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            pfx, lbound(pfx), ubound(pfx), &
            pfy, lbound(pfy), ubound(pfy), &
#if BL_SPACEDIM == 3
            pfz, lbound(pfz), ubound(pfz), &
#endif
            amrex_geom(lev)%dx, dt)
       
       if (do_reflux) then
          do idim = 1, amrex_spacedim
             pf => fluxes(idim)%dataptr(mfi)
             pfab => flux(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
             pf       (tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:) = &
                  pfab(tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:)
          end do
       end if

       ! advance particles on this tile
       particles => pc%get_particles(lev, mfi)
       call advect_particles(particles, size(particles), &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            dt, amrex_geom(lev)%dx, amrex_problo)
       
    end do
    call amrex_mfiter_destroy(mfi)
    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(uface(idim))
       call amrex_fab_destroy( flux(idim))
    end do
    !$omp end parallel

    if (do_reflux) then
       ! Note that the fluxes have already been scaled by dt and area.
       if (lev > 0) then
          call flux_reg(lev)%fineadd(fluxes, 1.0_amrex_real)
       end if

       if (lev < amrex_get_finest_level()) then
          call flux_reg(lev+1)%crseinit(fluxes, -1.0_amrex_real)
       end if

       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(fluxes(idim))
       end do
    end if

    call amrex_multifab_destroy(phiborder)
  end subroutine advance

end module evolve_module
