!> Miscellaneous functions for querying DVODE/CVODE integration statistics. Since DVODE and CVODE use different data structures to
!> track this data, we use our own data structure to capture the data from both solvers.

module integrator_stats_mod
  use, intrinsic :: iso_c_binding
  implicit none

  interface get_integrator_stats
    module procedure get_integrator_stats_dvode
    module procedure get_integrator_stats_cvode
  end interface get_integrator_stats

  type integrator_stats_t
    integer(c_long) :: n_steps
    integer(c_long) :: n_rhs_eval
    integer(c_int) :: last_order
    integer(c_int) :: current_order
    real(c_double) :: last_step_size
    real(c_double) :: current_step_size
    real(c_double) :: current_time
  end type integrator_stats_t

  contains

    function get_integrator_stats_cvode (cvmem) result (stats)
      use, intrinsic :: iso_c_binding
      use cvode_interface
      implicit none
      type(c_ptr), value :: cvmem
      type(integrator_stats_t) :: stats
      integer(c_int) :: ierr

      ierr = FCVodeGetNumSteps(cvmem, stats%n_steps)
      ierr = FCVodeGetNumRhsEvals(cvmem, stats%n_rhs_eval)
      ierr = FCVodeGetLastOrder(cvmem, stats%last_order)
      ierr = FCVodeGetCurrentOrder(cvmem, stats%current_order)
      ierr = FCVodeGetLastStep(cvmem, stats%last_step_size)
      ierr = FCVodeGetCurrentStep(cvmem, stats%current_step_size)
      ierr = FCVodeGetCurrentTime(cvmem, stats%current_time)

    end function get_integrator_stats_cvode

    pure function get_integrator_stats_dvode (iwork, rwork) result (stats)
      implicit none
      integer, intent(in) :: iwork(:)
      double precision, intent(in) :: rwork(:)
      type(integrator_stats_t) :: stats

      stats%n_steps = iwork(11)
      stats%n_rhs_eval = iwork(12)
      stats%last_order = iwork(14)
      stats%current_order = iwork(15)
      stats%last_step_size = rwork(11)
      stats%current_step_size = rwork(12)
      stats%current_time = rwork(13)

    end function get_integrator_stats_dvode

    subroutine print_integrator_stats(stats, print_header)
      implicit none

      type(integrator_stats_t), intent(in) :: stats
      logical, optional :: print_header

      if (present(print_header)) then
        if(print_header) then
          write(*, '(7a20)') 'n_steps', 'n_rhs_eval', 'last_order', &
            'current_order', 'last_step_size', 'current_step_size', &
            'current_time'
        end if
      end if

      write(*, '(4i20, 3es20.3e2)') stats%n_steps, stats%n_rhs_eval, &
        stats%last_order, stats%current_order, stats%last_step_size, &
        stats%current_step_size, stats%current_time

    end subroutine print_integrator_stats

end module integrator_stats_mod
