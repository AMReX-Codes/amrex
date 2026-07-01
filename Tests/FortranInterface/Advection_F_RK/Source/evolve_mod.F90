#include <AMReX_Config.H>

module evolve_module

  use iso_c_binding
  use amrex_amr_module

  implicit none
  private

  public :: evolve

  integer, parameter :: state_ngrow = 3
  type, bind(c) :: rk_context
     integer(c_int) :: lev, substep, ncycle, rk_order
     real(amrex_real) :: time, dt
  end type rk_context

contains

  subroutine evolve ()
    use my_amr_module, only : stepno, max_step, stop_time, dt, plot_int
    use amr_data_module, only : t_new
    use compute_dt_module, only : compute_dt
    use plotfile_module, only : writeplotfile
    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, substep, finest_level

    cur_time = t_new(0)
    last_plot_file_step = 0

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
    use amr_data_module, only : t_old, t_new, phi_old, phi_new, flux_reg
    use averagedown_module, only : averagedownto
    use fillpatch_module, only : reset_fillpatcher
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time

    integer, allocatable, save :: last_regrid_step(:)
    integer :: k, old_finest_level, finest_level, fine_substep

    if (regrid_int .gt. 0) then
       if (.not.allocated(last_regrid_step)) then
          allocate(last_regrid_step(0:amrex_max_level))
          last_regrid_step = 0
       end if

       if (lev .lt. amrex_max_level .and. stepno(lev) .gt. last_regrid_step(lev)) then
          if (mod(stepno(lev), regrid_int) .eq. 0) then
             old_finest_level = amrex_get_finest_level()
             call amrex_regrid(lev, time)
             finest_level = amrex_get_finest_level()

             do k = lev, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
             end do
          end if
       end if
    end if

    stepno(lev) = stepno(lev)+1

    t_old(lev) = time
    t_new(lev) = time + dt(lev)
    call amrex_multifab_swap(phi_old(lev), phi_new(lev))

    if (do_reflux .and. lev < amrex_get_finest_level()) then
       call flux_reg(lev+1)%setval(0.0_amrex_real)
    end if

    call advance(lev, time, dt(lev), stepno(lev), substep, nsubsteps(lev))

    if (lev .lt. amrex_get_finest_level()) then
       do fine_substep = 1, nsubsteps(lev+1)
          call timestep(lev+1, time+(fine_substep-1)*dt(lev+1), fine_substep)
       end do

       if (do_reflux) then
          call flux_reg(lev+1)%reflux(phi_new(lev), 1.0_amrex_real)
       end if

       call averagedownto(lev)
       call reset_fillpatcher(lev+1)
    end if

  end subroutine timestep

  subroutine advance (lev, time, dt, step, substep, ncycle)
    use my_amr_module, only : verbose, time_integrator
    integer, intent(in) :: lev, step, substep, ncycle
    real(amrex_real), intent(in) :: time, dt

    if (verbose .gt. 0 .and. amrex_parallel_ioprocessor()) then
       write(*,'(A, 1X, I0, 1X, A, 1X, I0, A, 1X, A, 1X, A, 1X, G0)') &
            "[Level", lev, "step", step, "] ADVANCE", time_integrator, "dt =", dt
    end if

    select case (time_integrator)
    case ("rk2")
       call advance_rk(lev, time, dt, substep, ncycle, 2)
    case ("rk3")
       call advance_rk(lev, time, dt, substep, ncycle, 3)
    case ("rk4")
       call advance_rk(lev, time, dt, substep, ncycle, 4)
    case default
       call amrex_abort("unknown myamr.time_integrator")
    end select
  end subroutine advance

  subroutine advance_rk (lev, time, dt, substep, ncycle, rk_order)
    use amrex_rungekutta_module, only : amrex_rungekutta_rk2, amrex_rungekutta_rk3, &
         amrex_rungekutta_rk4
    use amr_data_module, only : phi_old, phi_new
    integer, intent(in) :: lev, substep, ncycle, rk_order
    real(amrex_real), intent(in) :: time, dt

    type(rk_context), target :: ctx

    ctx%lev = lev
    ctx%substep = substep
    ctx%ncycle = ncycle
    ctx%rk_order = rk_order
    ctx%time = time
    ctx%dt = dt

    select case (rk_order)
    case (2)
       call amrex_rungekutta_rk2(phi_old(lev), phi_new(lev), time, dt, c_loc(ctx), &
            rk_rhs_callback, rk_fill_callback)
    case (3)
       call amrex_rungekutta_rk3(phi_old(lev), phi_new(lev), time, dt, c_loc(ctx), &
            rk_rhs_callback, rk_fill_callback, rk_store3_callback)
    case (4)
       call amrex_rungekutta_rk4(phi_old(lev), phi_new(lev), time, dt, c_loc(ctx), &
            rk_rhs_callback, rk_fill_callback, rk_store4_callback)
    case default
       call amrex_abort("unsupported RK order")
    end select
  end subroutine advance_rk

  subroutine rk_rhs_callback (ctx_ptr, stage, prhs, pstate, time, dtsub) &
       bind(c, name="amrex_advection_rk_rhs")
    type(c_ptr), value :: ctx_ptr, prhs, pstate
    integer(c_int), value :: stage
    real(amrex_real), value :: time, dtsub

    type(rk_context), pointer :: ctx
    type(amrex_multifab) :: rhs, state

    call c_f_pointer(ctx_ptr, ctx)
    rhs = prhs
    state = pstate

    call compute_rhs(int(ctx%lev), int(stage), state, rhs, time, dtsub, &
         int(ctx%ncycle), int(ctx%rk_order))
  end subroutine rk_rhs_callback

  subroutine rk_fill_callback (ctx_ptr, stage, pstate, time) &
       bind(c, name="amrex_advection_rk_fill")
    use fillpatch_module, only : fillpatch_stage
    type(c_ptr), value :: ctx_ptr, pstate
    integer(c_int), value :: stage
    real(amrex_real), value :: time

    type(rk_context), pointer :: ctx
    type(amrex_multifab) :: state

    call c_f_pointer(ctx_ptr, ctx)
    state = pstate

    call fillpatch_stage(int(ctx%lev), time, state, int(ctx%rk_order), int(stage), &
         int(ctx%substep), int(ctx%ncycle))
  end subroutine rk_fill_callback

  subroutine rk_store3_callback (ctx_ptr, ps_old, prk) &
       bind(c, name="amrex_advection_rk_store3")
    use fillpatch_module, only : store_rk3_fillpatcher
    type(c_ptr), value :: ctx_ptr, ps_old
    type(c_ptr), intent(in) :: prk(*)

    integer :: i
    type(rk_context), pointer :: ctx
    type(amrex_multifab) :: s_old, rk(3)

    call c_f_pointer(ctx_ptr, ctx)
    s_old = ps_old
    do i = 1, 3
       rk(i) = prk(i)
    end do

    call store_rk3_fillpatcher(int(ctx%lev), ctx%time, ctx%dt, s_old, rk, state_ngrow)
  end subroutine rk_store3_callback

  subroutine rk_store4_callback (ctx_ptr, ps_old, prk) &
       bind(c, name="amrex_advection_rk_store4")
    use fillpatch_module, only : store_rk4_fillpatcher
    type(c_ptr), value :: ctx_ptr, ps_old
    type(c_ptr), intent(in) :: prk(*)

    integer :: i
    type(rk_context), pointer :: ctx
    type(amrex_multifab) :: s_old, rk(4)

    call c_f_pointer(ctx_ptr, ctx)
    s_old = ps_old
    do i = 1, 4
       rk(i) = prk(i)
    end do

    call store_rk4_fillpatcher(int(ctx%lev), ctx%time, ctx%dt, s_old, rk, state_ngrow)
  end subroutine rk_store4_callback

  subroutine compute_rhs (lev, stage, state, rhs, time, dtsub, ncycle, rk_order)
    use my_amr_module, only : do_reflux
    use amr_data_module, only : phi_new, flux_reg
    use face_velocity_module, only : get_face_velocity
    use advect_module, only : advection_dudt
    integer, intent(in) :: lev, stage, ncycle, rk_order
    real(amrex_real), intent(in) :: time, dtsub
    type(amrex_multifab), intent(inout) :: state
    type(amrex_multifab), intent(inout) :: rhs

    integer :: ncomp, idim
    logical :: nodal(3)
    real(amrex_real) :: face_area
    type(amrex_multifab) :: fluxes(amrex_spacedim)
    type(amrex_mfiter) :: mfi
    type(amrex_box) :: bx, tbx
    real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pin,pout,pux,puy,puz,pfx,pfy,pfz, &
         pf, pfab
    type(amrex_fab) :: uface(amrex_spacedim)
    type(amrex_fab) :: flux(amrex_spacedim)

    ncomp = phi_new(lev)%ncomp()

    if (do_reflux) then
       do idim = 1, amrex_spacedim
          nodal = .false.
          nodal(idim) = .true.
          call amrex_multifab_build(fluxes(idim), phi_new(lev)%ba, phi_new(lev)%dm, ncomp, 0, nodal)
       end do
    end if

    !$omp parallel private(idim,mfi,bx,tbx,pin,pout,pux,puy,puz,pfx,pfy,pfz,pf,pfab,uface,flux)
    do idim = 1, amrex_spacedim
       call uface(idim)%reset_omp_private()
       call flux(idim)%reset_omp_private()
    end do
    call amrex_mfiter_build(mfi, phi_new(lev), tiling=.true.)
    do while(mfi%next())
       bx = mfi%tilebox()

       pin  => state%dataptr(mfi)
       pout => rhs%dataptr(mfi)

       do idim = 1, amrex_spacedim
          tbx = bx
          call tbx%nodalize(idim)
          call flux(idim)%resize(tbx,ncomp)
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

       call get_face_velocity(time, &
            pux, lbound(pux), ubound(pux), &
            puy, lbound(puy), ubound(puy), &
#if BL_SPACEDIM == 3
            puz, lbound(puz), ubound(puz), &
#endif
            amrex_geom(lev)%dx, amrex_problo)

       call advection_dudt(bx%lo, bx%hi, &
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
            amrex_geom(lev)%dx)

       if (do_reflux) then
          do idim = 1, amrex_spacedim
             pf => fluxes(idim)%dataptr(mfi)
             pfab => flux(idim)%dataptr()
             tbx = mfi%nodaltilebox(idim)
             pf       (tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:) = &
                  pfab(tbx%lo(1):tbx%hi(1),tbx%lo(2):tbx%hi(2),tbx%lo(3):tbx%hi(3),:)
          end do
       end if
    end do
    call amrex_mfiter_destroy(mfi)
    do idim = 1, amrex_spacedim
       call amrex_fab_destroy(uface(idim))
       call amrex_fab_destroy( flux(idim))
    end do
    !$omp end parallel

    if (do_reflux) then
       ! advection_dudt returns unscaled dudt and physical face fluxes.
       ! Match the C++ reflux pattern by passing the RK-stage dtsub*dA scale
       ! directly to the per-direction flux-register update.
       do idim = 1, amrex_spacedim
#if BL_SPACEDIM == 2
          face_area = amrex_geom(lev)%dx(3-idim)
#else
          if (idim == 1) then
             face_area = amrex_geom(lev)%dx(2)*amrex_geom(lev)%dx(3)
          else if (idim == 2) then
             face_area = amrex_geom(lev)%dx(1)*amrex_geom(lev)%dx(3)
          else
             face_area = amrex_geom(lev)%dx(1)*amrex_geom(lev)%dx(2)
          end if
#endif

          ! fineadd contributes this level as the fine side of flux_reg(lev),
          ! and crseinit contributes this level as the coarse side of
          ! flux_reg(lev+1). The registers are refluxed after the fine level
          ! has caught up to the coarse level.
          if (lev > 0) then
             call flux_reg(lev)%fineadd(fluxes(idim), idim, dtsub*face_area)
          end if

          if (lev < amrex_get_finest_level()) then
             call flux_reg(lev+1)%crseinit(fluxes(idim), idim, -dtsub*face_area, add=.true.)
          end if
       end do

       do idim = 1, amrex_spacedim
          call amrex_multifab_destroy(fluxes(idim))
       end do
    end if

  end subroutine compute_rhs

end module evolve_module
