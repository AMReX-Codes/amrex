module bl_prof_module

  use bl_types

  implicit none

  interface
     subroutine cpu_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine cpu_second
     subroutine wall_second(s)
       use bl_types
       real(kind=dp_t) :: s
     end subroutine wall_second
  end interface

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

  type, private :: timer_rec
     real(dp_t) :: cum  = 0.0_dp_t
     real(dp_t) :: max  = -Huge(1.0_dp_t)
     real(dp_t) :: min  = +Huge(1.0_dp_t)
     real(dp_t) :: cld  = 0.0_dp_t
     integer :: cnt = 0
     integer :: cmx = -Huge(1)
     integer :: cmn = +Huge(1)
  end type timer_rec

  type, private :: activation_n
     integer :: reg = BL_PROF_NOT_REG
     type(timer_rec) :: rec
     type(timer_rec) :: rec_global
     real(dp_t) :: strt = 0.0_dp_t
     logical :: running = .false.
     logical :: fixed   = .false.
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
  private s_activation
  private greater_d
  private benchmark
  private p_value
  private p_stop
  private p_start
  private s_debug
  private s_cksum

  logical, private, save :: wall = .true.

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
    use parallel
    use bl_error_module
    logical, intent(in), optional :: on
    if ( present(on) ) call bl_prof_set_state(on)
    if ( associated(timers) .or. associated(the_stack) ) then
       call bl_error("BL_PROF_INITIALIZE: must be called once only")
    end if
!!    if ( parallel_q() ) wall = .false.
    !! Hand initialize the data structures
    allocate(timers(1))
    timers(1)%reg = 1
    timers(1)%name = "boxlib"
    allocate(the_stack(1))
    stk_p = 1
    allocate(the_call_tree%children(0))
    the_stack(1)%a_p => the_call_tree
    the_stack(1)%a_p%reg = 1
    call p_start(the_stack(1)%a_p)
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
    use bl_error_module
    !! will be used to tear down the execution stack
    integer i
    if ( stk_p /= 1 ) call bl_error("BL_PROF_FINALIZE: stk_p :", stk_p)
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
    use bl_error_module
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
             call bl_error("BL_PROF_TIMER_BUILD: name already registered: ", name)
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
    use bl_error_module
    type(bl_prof_timer), intent(inout), target :: bpt
    if ( .not. bpt_on ) return
    if ( stk_p < 1 ) call bl_error("BL_PROF_TIMER_DESTROY: stack underflow: ", stk_p)
    if ( the_stack(stk_p)%a_p%running ) then
       call bl_prof_timer_stop(bpt)
    end if
    stk_p = stk_p - 1
  end subroutine bl_prof_timer_destroy

  subroutine p_start(a)
    use bl_error_module
    type(activation_n), intent(inout) :: a
    if ( a%running ) call bl_error("P_START: should not be be running before start")
    if ( wall ) then
       call wall_second(a%strt)
    else
       call cpu_second(a%strt)
    end if
    a%running = .true.
  end subroutine p_start

  subroutine p_stop(a)
    use bl_error_module
    type(activation_n), intent(inout) :: a
    real(dp_t) :: time
    if ( .not. a%running ) call bl_error("P_STOP: should be be running before start")
    if ( wall ) then
       call wall_second(time)
    else
       call cpu_second(time)
    end if
    time = time - a%strt
    a%strt = 0.0_dp_t
    a%rec%cum  = a%rec%cum + time
    a%rec%max  = max(a%rec%max, time)
    a%rec%min  = min(a%rec%min, time)
    a%rec%cnt  = a%rec%cnt + 1
    a%running = .false.
  end subroutine p_stop

  subroutine bl_prof_timer_start(bpt)
    type(bl_prof_timer), intent(inout) :: bpt

    if ( .not. bpt_on ) return

    call p_start(the_stack(stk_p)%a_p)

  end subroutine bl_prof_timer_start

  subroutine bl_prof_timer_stop(bpt)
    type(bl_prof_timer), intent(inout) :: bpt

    if ( .not. bpt_on ) return
    call p_stop(the_stack(stk_p)%a_p)
    the_stack(stk_p)%a_p%running = .false.

  end subroutine bl_prof_timer_stop

  subroutine bl_prof_glean(fname, note)
    use parallel
    use bl_IO_module
    use sort_d_module
    use bl_error_module
    character(len=*), intent(in) :: fname
    character(len=*), intent(in), optional :: note
    integer :: un
    character(len=8) :: date
    character(len=10) :: time
    integer :: i, ii, j
    real(dp_t), allocatable :: sm(:,:)
    integer, allocatable :: ism(:)
    integer, allocatable :: cksums(:)
    integer              :: cksum(1)
    logical              :: ok

    if ( stk_p /= 1 ) call bl_error("BL_PROF_GLEAN: stk_p :", stk_p)

    call p_stop(the_stack(1)%a_p)

    if ( parallel_ioprocessor() ) then
       un = unit_new()
       open(unit = un, file = trim(fname), form = "formatted", access = "sequential", &
            status = "replace", &
            action = "write")
    end if

    allocate(sm(size(timers),0:4),ism(size(timers)))
    sm(:,0) = 0.0_dp_t; sm(:,1) = 0.0_dp_t; sm(:,2) = 0.0_dp_t; sm(:,3) = -Huge(sm); sm(:,4) = +Huge(sm)

    allocate(cksums(parallel_nprocs()))
    cksum(1) = 0
    call s_cksum(the_call_tree, cksum(1))
    call parallel_gather(cksum, cksums, 1)
    ok = .true.
    if ( parallel_ioprocessor() ) then
       do i = 2, size(cksums)
          if ( cksums(i-1) /= cksums(i) ) then
             ok = .false.
             exit
          end if
       end do
    end if
    call parallel_bcast(ok)

    if ( .not. ok) then
       !
       ! Attempt to print out something to help track down why trees are not the same.
       !
       call sort(sm(:,2), ism, greater_d)
       do j = 0, parallel_nprocs()-1
          if (j == parallel_myproc()) then
             print*, 'glean on CPU# ', j, 'size(ism):', size(ism)
             do i = 1, size(ism)
                ii = ism(i)
                print*, trim(timers(ii)%name)
             end do
          end if
          call parallel_barrier()
       end do
       call bl_error("*** bl_prof_glean: proc trees are NOT identical !!!")
    end if

    call s_activation(the_call_tree, sm, local = .false.)

    if ( parallel_ioprocessor() ) then
       !! Print the summary informationreg
       write(unit = un, fmt = '("* GLOBAL")')
       write(unit = un, fmt = '("  NPROCS = ", i5,/)') parallel_nprocs()
       call sort(sm(:,2), ism, greater_d)
       write(unit = un, fmt = '("REGION",TR40,"COUNT",TR8,"TOTAL",TR22,"SELF",TR23,"MAX",TR10,"MIN")')
       do i = 1, size(ism)
          ii = ism(i)
          write(unit = un, fmt = '(A40,1x,I10,F13.3,TR13,F13.3,TR13,2F13.3)') trim(timers(ii)%name), int(sm(ii,0)), sm(ii,1:)
       end do
    end if

    if ( parallel_ioprocessor() ) then
       write(unit = un, fmt = '()')
       write(unit = un, fmt = &
            '("REGION",TR40,"COUNT",TR8,"TOTAL",TR8,"CHILD",TR9,"SELF",TR10,"AVG",TR10,"MAX",TR10,"MIN")')
    end if

    call p_activation(the_call_tree, un, 0, local = .false.)

    if ( parallel_ioprocessor() ) then
       close(unit = un)
    end if

    if ( parallel_nprocs() > 1 ) then
       do i = 0, parallel_nprocs() - 1
          if ( parallel_myproc() == i ) then
             un = unit_new()
             open(unit = un, file = trim(fname), form = "formatted", access = "sequential", &
                  position = "append", status = "old", &
                  action = "write")
             if ( i == 0 ) then
                write(unit = un, fmt = '("* LOCAL ", i5 )') i
             end if
             write(unit = un, fmt = '(/,/, "** PROCESSOR ", i5 )') i
             sm(:,1) = 0.0_dp_t; sm(:,2) = 0.0_dp_t; sm(:,3) = -Huge(sm); sm(:,4) = +Huge(sm)
             call s_activation(the_call_tree, sm, local = .true.)
             !! Print the summary information
             call sort(sm(:,2), ism, greater_d)
             write(unit = un, fmt = '("REGION",TR40,"COUNT",TR8,"TOTAL",TR22,"SELF",TR23,"MAX",TR10,"MIN")')
             do j = 1, size(ism)
                ii = ism(j)
                write(unit = un, fmt = '(A40,1x,I10,F13.3,TR13,F13.3,TR13,2F13.3)') &
                     trim(timers(ii)%name), int(sm(ii,0)), sm(ii,1:)
             end do
             write(unit = un, fmt = '()')
             write(unit = un, fmt = &
                  '("REGION",TR40,"COUNT",TR8,"TOTAL",TR8,"CHILD",TR9,"SELF",TR10,"AVG",TR10,"MAX",TR10,"MIN")')
             call p_activation(the_call_tree, un, 0, local = .true.)
             close(unit = un)
          end if
          call parallel_barrier()
       end do
    end if
  end subroutine bl_prof_glean

  function p_value(a, local) result(r)
    use parallel
    type(activation_n), intent(inout) :: a
    type(timer_rec) :: r
    logical, intent(in) :: local
    real(kind=dp_t), parameter :: MIL_SEC = 1.0e3_dp_t
    integer :: i
    if ( local ) then
       a%rec%cld = 0.0_dp_t
       do i = 1, size(a%children)
          a%rec%cld = a%rec%cld + a%children(i)%rec%cum
       end do
       r = a%rec
    else
       if ( a%fixed ) then
          r = a%rec_global
          goto 999
       end if
       r%cnt = a%rec%cnt
       a%rec%cld = 0.0_dp_t
       do i = 1, size(a%children)
          a%rec%cld = a%rec%cld + a%children(i)%rec%cum
       end do
       call parallel_reduce(r%cld, a%rec%cld, MPI_MAX)
       call parallel_reduce(r%cum, a%rec%cum, MPI_MAX)
       call parallel_reduce(r%max, a%rec%max, MPI_MAX)
       call parallel_reduce(r%min, a%rec%min, MPI_MIN)
       call parallel_reduce(r%cmx, a%rec%cnt, MPI_MAX)
       call parallel_reduce(r%cmn, a%rec%cnt, MPI_MIN)
       a%rec_global = r
       a%fixed = .true.
    end if

999 continue

    r%cum = r%cum*MIL_SEC
    r%max = r%max*MIL_SEC
    r%min = r%min*MIL_SEC
    r%cld = r%cld*MIL_SEC

  end function p_value

  function greater_d(a,b) result(r)
    logical :: r
    real(kind=dp_t), intent(in) :: a, b
    r = a > b
  end function greater_d

  recursive subroutine s_debug(a)
    use parallel
    type(activation_n), intent(inout) :: a
    integer :: i
    write(unit = 100 + parallel_myproc(), fmt=*) timers(a%reg)%name
    do i = 1, size(a%children)
       call s_debug(a%children(i))
    end do
  end subroutine s_debug

  recursive subroutine s_cksum(a, cksum)
    type(activation_n), intent(inout) :: a
    integer, intent(inout) :: cksum
    integer :: i
    call hash(timers(a%reg)%name, cksum)
    do i = 1, size(a%children)
       call s_cksum(a%children(i), cksum)
    end do
    contains
      subroutine hash(s, cksum)
        integer, intent(inout) :: cksum
        character(len=*), intent(in) :: s
        integer :: i
        do i = 1, len_trim(s)
           cksum = mod(64*cksum + ichar(s(i:i)), 1234567)
        end do
      end subroutine hash
  end subroutine s_cksum

  recursive subroutine s_activation(a, sm, local)
    type(activation_n), intent(inout) :: a
    real(dp_t), intent(inout) :: sm(:,0:)
    logical, intent(in) :: local
    real(dp_t) :: cum, self, cum_children
    integer :: i
    type(timer_rec) :: trec
    trec = p_value(a, local)
    cum = trec%cum
    cum_children = trec%cld
    self = cum - cum_children
    sm(a%reg,0) = sm(a%reg,0) + trec%cnt
    sm(a%reg,1) = sm(a%reg,1) + cum
    sm(a%reg,2) = sm(a%reg,2) + self
    sm(a%reg,3) = max(sm(a%reg,3),trec%max)
    sm(a%reg,4) = min(sm(a%reg,4),trec%min)
    do i = 1, size(a%children)
       call s_activation(a%children(i), sm, local)
    end do
  end subroutine s_activation

  recursive subroutine p_activation(a, un, skip, local)
    use parallel
    use bl_IO_module
    use sort_d_module
    type(activation_n), intent(inout) :: a
    integer, intent(in) :: un, skip
    logical, intent(in) :: local
    integer :: i, ii
    character(len=40) :: nm
    real(dp_t) :: cum_chidren, self, avg, cum
    real(dp_t), allocatable :: ctm(:)
    integer, allocatable :: itm(:)
    type(timer_rec) :: trec

    nm(1:skip) = repeat(' ', skip)
    nm(skip+1:) = trim(timers(a%reg)%name)
    trec = p_value(a, local)
    cum = trec%cum
    cum_chidren = trec%cld
    self = cum - cum_chidren
    avg  = self/trec%cnt
    if ( local .or. parallel_ioprocessor() ) then
       write(unit = un, fmt = '(A,1x,i10,6(F13.3))') nm, &
            trec%cnt, cum, cum_chidren, self, avg, trec%max, trec%min
    end if
    allocate(ctm(size(a%children)),itm(size(a%children)))
    do i = 1, size(a%children)
       trec = p_value(a%children(i), local)
       ctm(i) = trec%cum
    end do
    call sort(ctm, itm, greater_d)
    do i = 1, size(a%children)
       ii = itm(i)
       call p_activation(a%children(ii), un, skip+1, local)
    end do
  end subroutine p_activation

  subroutine print_stack()
    ! does nothing -- this is needed by the backtrace version
  end subroutine print_stack


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
