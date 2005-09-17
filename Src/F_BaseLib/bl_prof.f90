module bl_prof_module

  use bl_error_module

  implicit none

  !! A registry of unique numbers corresponding to the names of the
  !! timers
  integer, parameter, private :: BL_PROF_MAX_NAME_LENGTH = 64
  integer, parameter, private :: BL_PROF_NOT_REG = 0
  type bl_prof_reg
     integer :: reg = BL_PROF_NOT_REG
     character(len=BL_PROF_MAX_NAME_LENGTH) :: name = ''
  end type bl_prof_reg

  !! This is the object that lives in user space
  type bl_prof_timer
     integer :: reg = BL_PROF_NOT_REG
  end type bl_prof_timer

  type bl_prof_timer_record
     integer :: reg = BL_PROF_NOT_REG
     type(bl_prof_timer_record), pointer :: next => Null()
  end type bl_prof_timer_record

  type bl_prof_stack_n
     type(bl_prof_timer_record), pointer :: rec_p => Null()
     type(bl_prof_stack_n), pointer :: next => Null()
  end type bl_prof_stack_n

  type bl_prof_stack
     type(bl_prof_stack_n), pointer :: head => Null()
  end type bl_prof_stack

  interface pop
     module procedure bl_prof_stack_pop
  end interface

  interface push
     module procedure bl_prof_stack_push
  end interface

  interface front
     module procedure bl_prof_stack_front
  end interface

  interface build
     module procedure bl_prof_timer_build
  end interface

  interface start
     module procedure bl_prof_timer_start
  end interface

  interface stop
     module procedure bl_prof_timer_stop
  end interface

  integer, parameter, private :: BL_PROF_MAX_TIMERS = 1024
  type(bl_prof_reg), private, save :: timers(BL_PROF_MAX_TIMERS)

  !! By default, we don't profile should be relatively low-overhead if
  !! used judiciously
  logical, private, save :: bpt_on = .false.

  type(bl_prof_stack_n), target, private, save :: base_node
  type(bl_prof_stack), private, save :: the_stack

contains

  subroutine bl_prof_set_state(on)
    logical, intent(in) :: on
    bpt_on = on
  end subroutine bl_prof_set_state

  function bl_prof_get_state() result(r)
    logical :: r
    r = bpt_on
  end function bl_prof_get_state

  subroutine bl_prof_initialize(on)
    logical, intent(in), optional :: on
    if ( present(on) ) call bl_prof_set_state(on)
    the_stack%head => base_node
  end subroutine bl_prof_initialize

  subroutine bl_prof_finalize
    !! will be used to tear down the execution stack
  end subroutine bl_prof_finalize

  subroutine bl_prof_stack_pop(bps)
    type(bl_prof_stack), intent(in) :: bps
  end subroutine bl_prof_stack_pop

  subroutine bl_prof_stack_push(bps)
    type(bl_prof_stack), intent(in) :: bps
  end subroutine bl_prof_stack_push

  function bl_prof_stack_front(bps) result(r)
    type(bl_prof_timer_record) :: r
    type(bl_prof_stack), intent(in) :: bps
  end function bl_prof_stack_front

  subroutine bl_prof_timer_init
  end subroutine bl_prof_timer_init

  subroutine bl_prof_timer_build(bpt, name)
    type(bl_prof_timer), intent(inout) :: bpt
    character(len=*), intent(in) :: name
    integer :: i

    if ( .not. bpt_on ) return
    !! If not registered, then register
    if ( bpt%reg == BL_PROF_NOT_REG ) then
       do i = 1, BL_PROF_MAX_TIMERS
          if ( timers(i)%reg == BL_PROF_NOT_REG ) exit
          if ( timers(i)%name == name ) then
             call bl_error("BL_PROF_TIMER_BUILD: name already registered", name)
          end if
       end do
       if ( i > BL_PROF_MAX_TIMERS ) then
          call bl_error("BL_PROF_TIMER_BUILD: out of timers")
       end if
       timers(i)%reg = i
       timers(i)%name = name
       bpt%reg = i
    end if

  end subroutine bl_prof_timer_build

  subroutine bl_prof_timer_start(bpt)
    type(bl_prof_timer), intent(inout) :: bpt
    type(bl_prof_timer_record) :: rec
    type(bl_prof_timer_record), pointer :: recp
    if ( .not. bpt_on ) return
    rec = front(the_stack)
    recp => rec%next
    do while ( associated(recp) ) 
       if ( recp%reg == bpt%reg ) then
          exit
       end if
       recp => recp%next
    end do
    print *, 'associated ', associated(recp)
    if ( .not. associated(recp) ) then
       allocate(recp)
       recp%next => rec%next
       rec%next => recp
    end if

  end subroutine bl_prof_timer_start

  subroutine bl_prof_timer_stop(bpt)
    type(bl_prof_timer), intent(inout) :: bpt
    if ( .not. bpt_on ) return
  end subroutine bl_prof_timer_stop

  subroutine bl_prof_glean(fname, note)
    use bl_IO_module
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: note
    integer :: un
    character(len=8) :: date
    character(len=10) :: time
    integer :: i, cnt
    type(bl_prof_timer_record) :: rec

    un = unit_new()
    open(unit = un, file = trim(fname), &
         form = "formatted", access = "sequential", &
         status = "replace", action = "write")

    !! Preamble
    call date_and_time(date = date, time = time)
    write(unit = un, &
         fmt = '("(* BL PROF results for ", A2,"/",A2,"/",A4, " ", A2,":",A2,":",A2)') &
         date(5:6), date(7:8), date(1:4), &
         time(1:2), time(3:4), time(5:6)
    if ( present(note) ) then
       write(unit = un, fmt = '("   Note: ", A)') note
    end if
    write(unit = un, fmt = '("*)")')

    !!
    write(unit = un, fmt = '("registeredTimers = {")')
    cnt = 0
    do i = 1, BL_PROF_MAX_TIMERS
       if ( timers(i)%reg == BL_PROF_NOT_REG ) exit
       cnt = cnt + 1
    end do
    do i = 1, cnt
       write(unit = un, fmt = '("{", i0, ", ", A, "}")', advance = 'no') &
            timers(i)%reg, trim(timers(i)%name)
       if ( i < cnt ) then
          write(unit = un, fmt = '(",")')
       else
          write(unit = un, fmt = '()')
       end if
    end do
    write(unit = un, fmt = '("};")')
    !! Now grovel through the stack
    rec = bl_prof_stack_front(the_stack)

    !! Done
    close(unit = un)
  end subroutine bl_prof_glean

end module bl_prof_module

subroutine t_bl_prof
  use bl_prof_module
  implicit none
  type(bl_prof_timer), save :: bpt

  call bl_prof_initialize(on = .true.)
  call build(bpt, "t_bl_prof")
  call start(bpt)
  call t()
  call t()
  call stop(bpt)
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize

contains

  subroutine t
    type(bl_prof_timer), save :: bpt
    call build(bpt, "t")
    call start(bpt)
    call stop(bpt)
  end subroutine t

  subroutine g
    type(bl_prof_timer), save :: bpt
    call build(bpt, "t")        ! (sic)
    call start(bpt)
    call stop(bpt)
  end subroutine g

end subroutine t_bl_prof
