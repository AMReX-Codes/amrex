subroutine t_bl_prof
  use parallel
  use bl_prof_module
  implicit none
  type(bl_prof_timer), save :: bpt
  real :: a(100), b(100)
  integer :: i

  call bl_prof_initialize(on = .true.)
  call build(bpt, "t_bl_prof")
  b = (/(i,i=1,100)/)
  call t()
  do i = 1, 100
     a = sin(b)
  end do
  call t1()
  call destroy(bpt)
  call bl_prof_glean("bl_prof_res")
  call bl_prof_finalize

contains

  subroutine t
    type(bl_prof_timer), save :: bpt
    integer i
    call build(bpt, "t", no_start = .true.)
    call start(bpt)
    do i = 1, 100
       a = sin(b)
       a = sin(b)
    end do
    call destroy(bpt)
  end subroutine t

  subroutine t1
    type(bl_prof_timer), save :: bpt
    integer i
    call build(bpt, "t1", no_start = .true.)
    call start(bpt)
    do i = 1, 100*parallel_myproc()
       a = sin(b)
       a = sin(b)
    end do
    call destroy(bpt)
  end subroutine t1

  subroutine g
    type(bl_prof_timer), save :: bpt
    call build(bpt, "t")        ! (sic)
    call start(bpt)
    call stop(bpt)
    call destroy(bpt)
  end subroutine g

end subroutine t_bl_prof
