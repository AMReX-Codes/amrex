!! A module for timers
module bl_timer_module

  use bl_types

  implicit none

  integer,         private, parameter :: TIMER_NAME_MAX_LEN = 20
  real(kind=dp_t), private, parameter :: zero = 0.0_dp_t
  real(kind=dp_t), private, parameter :: MIL_SEC = 1.0e3_dp_t

  type timer
     character(len=TIMER_NAME_MAX_LEN) :: name = ""
     real(kind=dp_t) :: strt = zero
     real(kind=dp_t) :: time = zero
     real(kind=dp_t) :: cum_time = zero
     integer(kind=ll_t) :: cnt = 0_ll_t
     logical :: running = .FALSE.
     logical :: wall = .TRUE.
  end type timer

  interface
     subroutine cpu_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second
     subroutine wall_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second
     subroutine cpu_second_tick(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second_tick
     subroutine wall_second_tick(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second_tick
  end interface

  private :: cpu_second, wall_second
 
contains

  subroutine timer_build(tm, str, wall)
    type(timer), intent(out) :: tm
    character(len=*), intent(in) :: str
    logical, intent(in), optional :: wall
    tm%name = str
    tm%time = zero
    tm%strt = zero
    tm%cum_time = zero
    tm%cnt = 0_ll_t
    tm%running = .FALSE.
    if ( present(wall) ) tm%wall = wall
  end subroutine timer_build

  subroutine timer_start(tm)
    use bl_error_module
    type(timer), intent(inout) :: tm

    if ( tm%running ) &
         call bl_error("TIMER_START: should not be be running before start")
    if ( tm%wall ) then
       call wall_second(tm%strt)
    else
       call cpu_time(tm%strt)
    end if
    tm%running = .TRUE.
    tm%time = zero
  end subroutine timer_start

  subroutine timer_stop(tm)
    use bl_error_module
    type(timer), intent(inout) :: tm

    if ( .not. tm%running ) &
         call bl_error("TIMER_STOP: should be running before stop")
    if ( tm%wall ) then
       call wall_second(tm%time)
    else
       call cpu_time(tm%time)
    end if
    tm%time = tm%time-tm%strt
    tm%strt = zero
    tm%running = .FALSE.
    tm%cum_time = tm%cum_time + tm%time
    tm%cnt = tm%cnt + 1_ll_t

  end subroutine timer_stop

  ! This is dangerous since it forces a reduction!
  function timer_value(tm, total, proc) result(r)
    use parallel
    use bl_error_module
    real(kind=dp_t) :: r
    type(timer), intent(in) :: tm
    logical, intent(in), optional :: total
    integer, intent(in), optional :: proc
    logical :: ltotal
    real(kind=dp_t) :: r1
    ltotal = .FALSE.; if ( present(total) ) ltotal = total
    if ( tm%running ) &
         call bl_error("TIMER_VALUE: should be stopped before getting value")
    if ( ltotal ) then
       r1 = tm%cum_time
    else
       r1 = tm%time
    end if
    call parallel_reduce(r, r1, MPI_MAX, proc = proc)
    r = r1*MIL_SEC
  end function timer_value

  subroutine timer_print(tm, str, unit, advance, total)
    use parallel
    use bl_IO_module
    use bl_string_module
    type(timer), intent(in) :: tm
    character(len=*), intent(in), optional :: str
    character(len=*), intent(in), optional :: advance
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: total
    integer :: un
    character(len=3) :: adv
    logical :: ltotal
    real(dp_t) :: r

    !! Timer allows 9.9999999e9 as largest number of milliseconds print
    !! this is 3 years.
    !!      3301.5550
    !! F15.4 ( I have used a few timers with resolutions as low as .95 us.
    !!> 1234567890
    !!            1
    !!<            2345
    
    r = timer_value(tm, total, proc = parallel_ioprocessornode())
    if ( parallel_IOProcessor() ) then
       un = unit_stdout(unit)
       adv = unit_advance(advance)
       ltotal = .false.; if ( present(total) ) ltotal = total
       if ( present(str) ) then
          write(unit=un, fmt='("(*", A, "*)")', advance = adv) str
       end if
       if ( ltotal ) then
          write(unit=un, &
               fmt='("timer[",A,",",i15,",",F15.4,"]")', advance = 'no') &
               adjustr(tm%name), tm%cnt, tm%cum_time*MIL_SEC
       else
          write(unit=un, fmt='("timer[",A,",",F15.4,"]")', advance = 'no') &
               adjustr(tm%name), tm%time*MIL_SEC
       end if
       if ( eq_i(adv, 'YES') ) write(unit = un, fmt='()')
    end if

  end subroutine timer_print

  function timer_tick() result(r)
    real(kind=dp_t) :: r
    call wall_second_tick(r)
  end function timer_tick

  function MY_CPU_SECOND_TICK() result(r)
    real(kind=dp_t) :: r
    real(kind=dp_t), save ::  tickval = -1
    real(kind=dp_t) :: t1, t2
    integer :: cnt, icnt
    if ( tickval < 0 ) then
       tickval = Huge(tickval)
       do icnt = 1, 1000
          call cpu_time(t1)
          do cnt = 1, 1000
             call cpu_time(t2)
             if (t2 > t1) then
                tickval = min(tickval,t2 - t1)
                exit
             end if
          end do
       end do
       tickval = tickval
    end if
    r = tickval
  end function MY_CPU_SECOND_TICK

  function MY_WALL_SECOND_TICK() result(r)
    real(kind=dp_t) :: r
    real(kind=dp_t), save ::  tickval = -1
    real(kind=dp_t) :: t1, t2
    integer :: cnt, icnt
    if ( tickval < 0 ) then
       tickval = Huge(tickval)
       do icnt = 1, 1000
          call wall_second(t1)
          do cnt = 1, 1000
             call wall_second(t2)
             if (t2 > t1) then
                tickval = min(tickval,t2 - t1)
                exit
             end if
          end do
       end do
       tickval = tickval
    end if
    r = tickval
  end function MY_WALL_SECOND_TICK

end module bl_timer_module
