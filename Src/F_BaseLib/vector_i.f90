module vector_i_module
  
  implicit none

  type vector_i
     private
     integer :: end  = 0
     integer :: size = 0
     integer :: bump = 2
     integer, dimension(:), pointer :: d => NULL()
  end type vector_i

  interface capacity
    module procedure vector_capacity_i
  end interface capacity

  interface size
     module procedure vector_size_i
  end interface size

  interface build
     module procedure vector_build_i
     module procedure vector_build_v_i
  end interface build

  interface destroy
     module procedure vector_destroy_i
  end interface destroy

  interface empty
     module procedure vector_empty_i
  end interface empty

  interface at
     module procedure vector_at_i
  end interface at
  interface set_at
     module procedure vector_at_set_i
  end interface set_at

  interface dataptr
    module procedure vector_dataptr_i
  end interface dataptr

  interface push_back
     module procedure vector_push_back_i
  end interface push_back

  interface pop_back
     module procedure vector_pop_back_i
  end interface pop_back

  interface front
     module procedure vector_front_i
  end interface front

  interface erase
     module procedure vector_erase_i
  end interface erase

  interface back
     module procedure vector_back_i
  end interface back

  interface reserve
    module procedure vector_reserve_i
  end interface reserve

  interface resize
     module procedure vector_resize_i
  end interface resize

  interface reverse
     module procedure vector_reverse_i
  end interface reverse

  interface swap
     module procedure vector_swap_i
  end interface swap

  interface clear
     module procedure vector_clear_i
  end interface clear

  interface insert
     module procedure vector_insert_i
     module procedure vector_insert_n_i
  end interface insert

  interface print
     module procedure vector_print_i
  end interface print

  interface operator (.EQ.)
     module procedure vector_equal_i
  end interface
  interface equal
     module procedure vector_equal_i
  end interface equal
  interface operator (.NE.)
     module procedure vector_not_equal_i
  end interface
  interface not_equal
     module procedure vector_not_equal_i
  end interface not_equal

contains

  subroutine vector_build_v_i(d, values)
    type(vector_i), intent(out) :: d
    integer, intent(in), dimension(:) :: values

    d%size = size(values)
    allocate(d%d(d%size))
    d%d = values
    d%end  = size(values)

  end subroutine vector_build_v_i

  subroutine vector_build_i(d, size, value, bump)
    type(vector_i), intent(out) :: d
    integer, intent(in), optional :: size
    integer, intent(in), optional :: value
    integer, intent(in), optional :: bump

    integer :: v
        v = 0
    if ( present(value) ) v = value
    if ( present(size) )  d%size = size
    if ( present(bump) )  d%bump = bump
    allocate(d%d(d%size))
    d%d(1:d%size) = v

  end subroutine vector_build_i

  subroutine vector_destroy_i(d)
    type(vector_i), intent(inout) :: d

    call clear(d)
    d%bump=2

  end subroutine vector_destroy_i

  function vector_empty_i(d) result(r)
    logical :: r
    type(vector_i), intent(in) :: d

    r = d%size == 0

  end function vector_empty_i

  function vector_size_i(d) result(r)
    integer :: r
    type(vector_i), intent(in) :: d

    r = d%size

  end function vector_size_i

  function vector_capacity_i(d) result(r)
    integer :: r
    type(vector_i), intent(in) :: d

    r = size(d%d)

  end function vector_capacity_i

  function vector_at_i(d, i) result(r)
    integer :: r
    integer, intent(in) :: i
    type(vector_i), intent(in) :: d

    r = d%d(i)

  end function vector_at_i

  subroutine vector_at_set_i(d, i, v)
    integer :: v
    integer, intent(in) :: i
    type(vector_i), intent(inout) :: d

    d%d(i) = v

  end subroutine vector_at_set_i

  function vector_dataptr_i(d,lo,hi) result(r)
    integer, pointer, dimension(:) :: r
    type(vector_i), intent(in) :: d
    integer, intent(in), optional :: lo
    integer, intent(in), optional :: hi

    if ( present(lo) .AND. present(hi) ) then
       r => d%d(lo:hi)
    else if ( present(lo) ) then
       r => d%d(lo:d%end)
    else if ( present(hi) ) then
       r => d%d(1:hi)
    else
       r => d%d(1:d%end)
    end if

  end function vector_dataptr_i

  subroutine vector_resize_i(d, size, value)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: size
    integer, intent(in), optional ::  value
    integer, dimension(:), pointer :: np
    integer :: v

        v = 0
    if ( present(value) ) v = value
    if ( size <= vector_capacity_i(d) ) then
       d%d(1:d%size) = d%d(1:d%end)
       d%d(d%size+1:size) = v
    else
       allocate(np(size))
       np(1:d%size) = d%d(1:d%end)
       np(d%size+1:size) = v
       if ( associated(d%d) ) then
          deallocate(d%d)
       end if
       d%d => np
    end if
    d%size  = size
    d%end   = size

  end subroutine vector_resize_i

  subroutine vector_set_bump_i(d, bump)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: bump

    d%bump = bump

  end subroutine vector_set_bump_i

  subroutine vector_reserve_i(d, size)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: size
    integer, dimension(:), pointer :: np

    if ( size <= vector_capacity_i(d) ) return
    allocate(np(size))
    np(1:d%size) = d%d(1:d%size)
    if ( associated(d%d) ) then
       deallocate(d%d)
    end if
    d%d => np

  end subroutine vector_reserve_i

  subroutine vector_shrink_wrap_i(d)
    type(vector_i), intent(inout) :: d
    integer, dimension(:), pointer :: np

    if ( size(d%d) == d%size ) return
    allocate(np(d%size))
    np(1:d%size) = d%d(1:d%size)
    if ( associated(d%d) ) then
       deallocate(d%d)
    end if
    d%d => np
    d%end   = d%size

  end subroutine vector_shrink_wrap_i

  subroutine vector_push_back_i(d,v)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: v

    if ( d%size >= vector_capacity_i(d) ) then
       call vector_reserve_i(d,max(d%size+1,d%size*d%bump))
    end if
    d%size = d%size + 1
    d%end = d%end + 1
    d%d(d%end) = v

  end subroutine vector_push_back_i

  subroutine vector_pop_back_i(d)
    type(vector_i), intent(inout) :: d

    d%size = d%size - 1
    d%end = d%end - 1

  end subroutine vector_pop_back_i

  function vector_front_i(d) result(r)
    type(vector_i), intent(in) :: d

    integer :: r
    r = d%d(1)

  end function vector_front_i

  function vector_back_i(d) result(r)
    type(vector_i), intent(in) :: d

    integer :: r
    r = d%d(d%end)

  end function vector_back_i

  subroutine vector_reverse_i(d)
    type(vector_i), intent(inout) :: d
    integer :: i,n

    n = size(d)
    do i = 1, n/2
       call sw(d%d(i),d%d(n-i+1))
    end do

  contains
    subroutine sw(a,b)
      integer, intent(inout) :: a, b
      integer :: t

      t = a; a = b; b = t

    end subroutine sw

  end subroutine vector_reverse_i

  subroutine vector_swap_i(a,b)
    type(vector_i), intent(inout) :: a,b
    type(vector_i) :: t

    t = a; a = b; b = t

  end subroutine vector_swap_i

  function vector_erase_i(d,i) result(r)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: i
    integer :: r

    d%d(i:d%size-1) = d%d(i+1:d%size)
    d%size = d%size - 1
    d%end =  d%end  - 1
    r = i
  end function vector_erase_i

  subroutine vector_insert_i(d, i, elem)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: i

    integer, intent(in) :: elem
    if ( d%size >= vector_capacity_i(d) ) then
       call vector_reserve_i(d,max(d%size+1,d%size*d%bump))
    end if
    d%size = d%size + 1
    d%end  = d%end + 1
    d%d(i+1:d%size) = d%d(i:d%size-1)
    d%d(i) = elem

  end subroutine vector_insert_i

  subroutine vector_insert_n_i(d, i, elem, n)
    type(vector_i), intent(inout) :: d
    integer, intent(in) :: n
    integer, intent(in) :: i
    integer, intent(in) :: elem
    integer j

    do j = 1, n
       call vector_insert_i(d,i,elem)
    end do

  end subroutine vector_insert_n_i

  subroutine vector_clear_i(d)
    type(vector_i), intent(inout) :: d

    d%size = 0
    d%end =  0
    if ( associated(d%d) ) then
       deallocate(d%d)
    end if

  end subroutine vector_clear_i

  function vector_equal_i(d1, d2) result(r)
    logical :: r
    type(vector_i), intent(in) :: d1, d2
    integer :: i

    r = .TRUE.
    if ( size(d1) .ne. size(d2) ) then
       r = .FALSE.
       return
    end if
    do i = 1, size(d1)
       if ( d1%d(i) .NE. d2%d(i) ) then
          r = .FALSE.
          exit
       end if
    end do
  end function vector_equal_i
  function vector_not_equal_i(d1, d2) result(r)
    logical :: r
    type(vector_i), intent(in) :: d1, d2
    r = .not. equal(d1,d2)
  end function vector_not_equal_i

  subroutine vector_print_i(d, str, unit)
    use bl_IO_module
    type(vector_i), intent(in) :: d
    character (len=*), intent(in), optional :: str
    integer, intent(in), optional :: unit
    integer :: un

    un = unit_stdout(unit)
    if ( present(str) ) write(unit=un, fmt='(A)', advance='no') str
    if ( empty(d) ) then
       write(unit=un, fmt='("Empty")')
    else
       write(unit=un, fmt=*) d%d(1:d%end)
    end if

  end subroutine vector_print_i

end module vector_i_module
