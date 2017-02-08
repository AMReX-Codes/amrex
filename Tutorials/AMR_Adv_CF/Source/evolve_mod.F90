module evolve_module

  use amrex_amr_module

  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()
    use my_amr_module, only : phi_old, phi_new, stepno, max_step, t_new, stop_time, dt, plot_int
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
    use my_amr_module, only : regrid_int, stepno, nsubsteps, last_regrid_step, dt, do_reflux, &
         t_old, t_new, phi_old, phi_new
    use averagedown_module, only : averagedownto
    integer, intent(in) :: lev, substep
    real(amrex_real), intent(in) :: time
    
    integer :: i, k, old_finest_level, finest_level, fine_substep

    if (regrid_int .ge. 0) then
       do i = lev, amrex_max_level-1
          finest_level = amrex_get_finest_level()
          if (i .gt. finest_level) exit

          if (stepno(i) .gt. last_regrid_step(i) .and. mod(stepno(i), regrid_int) .eq. 0) then
             old_finest_level = finest_level
             call amrex_regrid(i, time)
             ! note that finest level can change during regrid
             finest_level = amrex_get_finest_level()

             do k = i, finest_level
                last_regrid_step(k) = stepno(k)
             end do

             do k = old_finest_level+1, finest_level
                dt(k) = dt(k-1) / amrex_ref_ratio(k-1)
             end do
          end if
       end do
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
          if (amrex_parallel_ioprocessor()) then
             print *, "TODO: reflux"
          end if
       end if

       call averagedownto(lev)
    end if    
  end subroutine timestep

  ! Given phi_old(lev), compute phi_new(lev)
  subroutine advance (lev, time, dt, step, substep, nsub)
    use my_amr_module, only : verbose, phi_new, phi_old
    integer, intent(in) :: lev, step, substep, nsub
    real(amrex_real), intent(in) :: time, dt

    integer, parameter :: ngrow = 3
    integer :: ncomp
    real(amrex_real) :: ctr_time
    type(amrex_multifab) :: phiborder

    if (verbose .gt. 0 .and. amrex_parallel_ioprocessor()) then
       write(*,'(A, 1X, I0, 1X, A, 1X, I0, A, 1X, G0)') &
            "[Level", lev, "step", step, "] ADVANCE with dt =", dt
    end if

    ncomp = phi_new(lev)%ncomp()

    call amrex_multifab_build(phiborder, phi_new(lev)%ba, phi_new(lev)%dm, ncomp, ngrow)

!    call fillpatch(lev, time, 

    ctr_time = time + 0.5_amrex_real * dt

    call amrex_multifab_destroy(phiborder)

  end subroutine advance

end module evolve_module
