module evolve_module

  use amrex_amr_module

  implicit none
  private

  public :: evolve

contains

  subroutine evolve ()
    use my_amr_module, only : phi_old, phi_new, istep, max_step, t_new, stop_time, dt, plot_int
    use compute_dt_module, only : compute_dt
    use plotfile_module, only : writeplotfile
    real(amrex_real) :: cur_time
    integer :: last_plot_file_step, step, lev, iteration, finest_level

    cur_time = t_new(0)
    last_plot_file_step = 0;
    
    do step = istep(0), max_step-1
       if (cur_time .ge. stop_time) exit

       if (amrex_parallel_ioprocessor()) then
          print *, ""
          print *, "STEP", step+1, "starts ..."
       end if

       call compute_dt()

       lev = 0
       iteration = 1
!       call timestep(lev, cur_time, iteration)
 
       cur_time = cur_time + dt(0)

       if (amrex_parallel_ioprocessor()) then
          print *, "STEP", step+1, "end. TIME =", cur_time, "DT =", dt(0)
       end if

       ! sync up time
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

    if (plot_int .gt. 0 .and. istep(0) .gt. last_plot_file_step) then
       call writeplotfile()
    end if

  end subroutine evolve

end module evolve_module
