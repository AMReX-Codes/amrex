module bl_prof_module

  implicit none

  !! This is the object that lives in user space
  type bl_prof_timer
     integer :: reg = -1
  end type bl_prof_timer

  interface build
     module procedure bl_prof_timer_build
  end interface

  interface stop
     module procedure bl_prof_timer_stop
  end interface

  interface start
     module procedure bl_prof_timer_start
  end interface

  interface destroy
     module procedure bl_prof_timer_destroy
  end interface

contains

  subroutine bl_prof_set_state(on)
    logical, intent(in) :: on
  end subroutine bl_prof_set_state

  function bl_prof_get_state() result(r)
    logical :: r
    r = .false.
  end function bl_prof_get_state

  subroutine bl_prof_initialize(on)
    logical, intent(in), optional :: on
  end subroutine bl_prof_initialize

  subroutine bl_prof_finalize
  end subroutine bl_prof_finalize

  subroutine bl_prof_timer_init
  end subroutine bl_prof_timer_init

  subroutine bl_prof_timer_build(bpt, name, no_start)
    type(bl_prof_timer), intent(inout), target :: bpt
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: no_start
  end subroutine bl_prof_timer_build

  subroutine bl_prof_timer_destroy(bpt)
    type(bl_prof_timer), intent(inout), target :: bpt
  end subroutine bl_prof_timer_destroy

  subroutine bl_prof_timer_start(bpt)
    type(bl_prof_timer), intent(inout) :: bpt
  end subroutine bl_prof_timer_start

  subroutine bl_prof_timer_stop(bpt)
    type(bl_prof_timer), intent(inout) :: bpt
  end subroutine bl_prof_timer_stop

  subroutine bl_prof_glean(fname, note)
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: note
  end subroutine bl_prof_glean

  subroutine print_stack()
    ! does nothing -- this is needed by the backtrace version
  end subroutine print_stack

end module bl_prof_module
