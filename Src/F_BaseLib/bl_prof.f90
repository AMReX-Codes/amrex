module bl_prof_module

  use bl_error_module
  use bl_timer_module

  implicit none

  !! A registry of unique numbers corresponding to the names of the
  !! timers
  integer, parameter, private :: BL_PROF_MAX_NAME_LENGTH = 64
  integer, parameter, private :: BL_PROF_NOT_REG = 0
  type, private :: bl_prof_reg
     integer :: reg = BL_PROF_NOT_REG
     character(len=BL_PROF_MAX_NAME_LENGTH) :: name = ''
     type(bl_prof_timer), pointer :: bpt => Null()
     logical :: managed = .false.
  end type bl_prof_reg

  !! This is the object that lives in user space
  type bl_prof_timer
     integer :: reg = BL_PROF_NOT_REG
  end type bl_prof_timer

  type, private :: activation_n
     integer :: reg = BL_PROF_NOT_REG
     type(timer) :: rec
     logical :: running = .false.
     type(activation_n), pointer :: children(:) => Null()
  end type activation_n

  type, private :: stack_n
     type(activation_n), pointer :: a_p
  end type stack_n

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

  type(bl_prof_reg), private, save, pointer :: timers(:)

  !! By default, we don't profile should be relatively low-overhead if
  !! used judiciously
  logical, private, save :: bpt_on = .false.

  type(stack_n), private, save, pointer :: the_stack(:) => Null()
  integer, private :: stk_p = 0

  type(activation_n), private, save, target :: the_call_tree

  private f_activation
  private p_activation
  private t_activation
  private s_activation
  private greater_d
  private benchmark

contains

  subroutine bl_prof_set_state(on)
    logical, intent(in) :: on
    bpt_on = on
  end subroutine bl_prof_set_state

  function bl_prof_get_state() result(r)
    logical :: r
    r = bpt_on
  end function bl_prof_get_state

  subroutine benchmark
    type(bl_prof_timer), save :: bpt
    call build(bpt, "bl_prof_benchmark")
    call destroy(bpt)
  end subroutine benchmark

  subroutine bl_prof_initialize(on)
    logical, intent(in), optional :: on
    if ( present(on) ) call bl_prof_set_state(on)
    if ( associated(timers) .or. associated(the_stack) ) then
       call bl_error("BL_PROF_INITIALIZE: must be called once only")
    end if
    !! Hand initialize the data structures
    allocate(timers(1))
    timers(1)%reg = 1
    timers(1)%name = "boxlib"
    allocate(the_stack(1))
    stk_p = 1
    allocate(the_call_tree%children(0))
    the_stack(1)%a_p => the_call_tree
    the_stack(1)%a_p%reg = 1
    call timer_start(the_stack(1)%a_p%rec)
    the_stack(1)%a_p%running = .true.
    call benchmark()
  end subroutine bl_prof_initialize

  recursive subroutine f_activation(a)
    type(activation_n), intent(inout) :: a
    integer i
    do i = 1, size(a%children)
       call f_activation(a%children(i))
    end do
    deallocate(a%children)
  end subroutine f_activation

  subroutine bl_prof_finalize
    !! will be used to tear down the execution stack
    integer i
    if ( stk_p /= 1 ) then
       call bl_error("BL_PROF_FINALIZE: stk_p :", stk_p)
    end if
    do i = 1, size(timers)
       if ( timers(i)%managed ) then
          deallocate(timers(i)%bpt)
       end if
    end do
    deallocate(timers)
    call f_activation(the_call_tree)
    deallocate(the_stack)
  end subroutine bl_prof_finalize

  subroutine bl_prof_timer_init
  end subroutine bl_prof_timer_init

  subroutine bl_prof_timer_build(bpt, name, no_start)
    type(bl_prof_timer), intent(inout), target :: bpt
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: no_start
    integer :: i
    type(bl_prof_reg), pointer :: new_timers(:)
    type(stack_n), pointer :: new_stack(:)
    type(activation_n), pointer :: new_children(:), a_p
    logical :: lstart

    if ( .not. bpt_on ) return
    !! reverse the sense
    lstart = .true.; if ( present(no_start) ) lstart = .not. no_start
    !! If not registered, then register
    if ( bpt%reg == BL_PROF_NOT_REG ) then
       do i = 1, size(timers)
          if ( timers(i)%reg == BL_PROF_NOT_REG ) exit
          if ( timers(i)%name == name ) then
             call bl_error("BL_PROF_TIMER_BUILD: name already registered", name)
          end if
       end do
       if ( i > size(timers) ) then
          allocate(new_timers(i))
          new_timers(1:i-1) = timers(1:i-1)
          deallocate(timers)
          timers => new_timers
       end if
       timers(i)%reg = i
       timers(i)%name = name
       timers(i)%bpt => bpt
       bpt%reg = i
    end if

    !! Call tree placement
    do i = 1, size(the_stack(stk_p)%a_p%children)
       a_p => the_stack(stk_p)%a_p%children(i)
       if ( a_p%reg == bpt%reg ) exit
    end do
    if ( i > size(the_stack(stk_p)%a_p%children) ) then
       allocate(new_children(i))
       new_children(1:i-1) = the_stack(stk_p)%a_p%children(1:i-1)
       deallocate(the_stack(stk_p)%a_p%children)
       the_stack(stk_p)%a_p%children => new_children
       a_p => the_stack(stk_p)%a_p%children(i)
       allocate(a_p%children(0))
       a_p%reg = bpt%reg
    end if

    !! Stack adjustment
    if ( stk_p == size(the_stack) ) then
       allocate(new_stack(stk_p+1))
       new_stack(1:stk_p) = the_stack(1:stk_p)
       deallocate(the_stack)
       the_stack => new_stack
    else if ( stk_p > size(the_stack) ) then
       call bl_error("BL_PROF_TIMER_BUILD: stack super-overflow", stk_p-size(the_stack))
    end if
    stk_p = stk_p + 1
    the_stack(stk_p)%a_p => a_p

    if ( lstart ) then
       call bl_prof_timer_start(bpt)
    end if

  end subroutine bl_prof_timer_build

  subroutine bl_prof_timer_destroy(bpt)
    type(bl_prof_timer), intent(inout), target :: bpt
    if ( .not. bpt_on ) return
    if ( stk_p < 1 ) then
       call bl_error("BL_PROF_TIMER_DESTROY: stack underflow: ", stk_p)
    end if
    if ( the_stack(stk_p)%a_p%running ) then
       call bl_prof_timer_stop(bpt)
    end if
    stk_p = stk_p - 1
  end subroutine bl_prof_timer_destroy

  subroutine bl_prof_timer_start(bpt)
    type(bl_prof_timer), intent(inout) :: bpt

    if ( .not. bpt_on ) return
    the_stack(stk_p)%a_p%running = .true.
    call timer_start(the_stack(stk_p)%a_p%rec)

  end subroutine bl_prof_timer_start

  subroutine bl_prof_timer_stop(bpt)
    type(bl_prof_timer), intent(inout) :: bpt

    if ( .not. bpt_on ) return
    call timer_stop(the_stack(stk_p)%a_p%rec)
    the_stack(stk_p)%a_p%running = .false.

  end subroutine bl_prof_timer_stop

  subroutine bl_prof_glean(fname, note)
    use bl_IO_module
    use sort_d_module
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: note
    integer :: un
    character(len=8) :: date
    character(len=10) :: time
    integer :: i, ii
    real(dp_t), allocatable :: sm(:,:)
    integer, allocatable :: ism(:)

    if ( stk_p /= 1 ) then
       call bl_error("BL_PROF_GLEAN: stk_p :", stk_p)
    end if
    call timer_stop(the_stack(1)%a_p%rec)

    if ( parallel_ioprocessor() ) then
       un = unit_new()
       open(unit = un, file = trim(fname), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
    end if

    allocate(sm(size(timers),2),ism(size(timers)))
    sm = 0.0_dp_t
    call s_activation(the_call_tree, sm)

    if ( parallel_ioprocessor() ) then
       !! Print the summary information
       call sort(sm(:,2), ism, greater_d)
       write(unit = un, fmt = '("REGION",TR24,"TOTAL", TR11, "SELF")')
       do i = 1, size(ism)
          ii = ism(i)
          write(unit = un, fmt = '(a20,2F15.3)') trim(timers(ii)%name), sm(ii,:)
       end do
    end if

    if ( parallel_ioprocessor() ) then
       write(unit = un, fmt = &
            '("REGION",TR20,"COUNT",TR10,"TOTAL", TR10, "CHILD", TR11, "SELF", TR12, "AVG")')
    end if

    call p_activation(the_call_tree, un, 0)
    if ( parallel_ioprocessor() ) then
       close(unit = un)
    end if
  end subroutine bl_prof_glean

  function t_activation(a) result(r)
    real(dp_t) :: r
    type(activation_n), intent(in) :: a
    integer :: i
    r = 0.0_dp_t
    do i = 1, size(a%children)
       r = r + timer_value(a%children(i)%rec, total = .true.)
    end do
  end function t_activation

  function greater_d(a,b) result(r)
    logical :: r
    real(kind=dp_t), intent(in) :: a, b
    r = a > b
  end function greater_d

  recursive subroutine s_activation(a, sm)
    type(activation_n), intent(in) :: a
    real(dp_t), intent(inout) :: sm(:,:)
    real(dp_t) :: cum, self, cum_children
    integer :: i
    cum = timer_value(a%rec, total = .true.)
    cum_children = t_activation(a)
    self = cum - cum_children
    sm(a%reg,1) = sm(a%reg,1) + cum
    sm(a%reg,2) = sm(a%reg,2) + self
    do i = 1, size(a%children)
       call s_activation(a%children(i), sm)
    end do
  end subroutine s_activation

  recursive subroutine p_activation(a, un, skip)
    use bl_IO_module
    use sort_d_module
    type(activation_n), intent(in) :: a
    integer, intent(in) :: un, skip
    integer :: i, cnt, ii
    character(len=20) :: nm
    real(dp_t) :: cum, cum_chidren, self, avg
    real(dp_t), allocatable :: ctm(:)
    integer, allocatable :: itm(:)
    
    nm(1:skip) = repeat(' ', skip)
    nm(skip+1:) = trim(timers(a%reg)%name)
    cnt = a%rec%cnt
    cum = timer_value(a%rec, total = .true.)
    cum_chidren = t_activation(a)
    self = cum - cum_chidren
    avg  = self/cnt
    if ( parallel_ioprocessor() ) then
       write(unit = un, fmt = '(A,1x,i10,4(F15.3))') nm, &
            cnt, cum, cum_chidren, self, avg
    end if
    allocate(ctm(size(a%children)),itm(size(a%children)))
    do i = 1, size(a%children)
       ctm(i) = timer_value(a%children(i)%rec, total = .true.)
    end do
    call sort(ctm, itm, greater_d)
    do i = 1, size(a%children)
       ii = itm(i)
       call p_activation(a%children(ii), un, skip+1)
    end do
  end subroutine p_activation

end module bl_prof_module

!! The C interface, not working
subroutine bl_prof_build_c(reg, name, n)
  use bl_prof_module
  use bl_string_module
  implicit none
  integer, intent(out) :: reg
  integer, intent(in) :: n
  integer, intent(in) :: name(n)
  character(len=n) :: cname
  call int2str(cname, name, n)
  reg = 0
end subroutine bl_prof_build_c

subroutine t_bl_prof
  use bl_prof_module
  implicit none
  type(bl_prof_timer), save :: bpt

  call bl_prof_initialize(on = .true.)
  call build(bpt, "t_bl_prof")
  call t()
  call t()
  call destroy(bpt)
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize

contains

  subroutine t
    type(bl_prof_timer), save :: bpt
    call build(bpt, "t", no_start = .true.)
    call start(bpt)
    call destroy(bpt)
  end subroutine t

  subroutine g
    type(bl_prof_timer), save :: bpt
    call build(bpt, "t")        ! (sic)
    call start(bpt)
    call stop(bpt)
    call destroy(bpt)
  end subroutine g

end subroutine t_bl_prof
