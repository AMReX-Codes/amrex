module vector_i_module
  
  implicit none

  real, private, parameter :: BUMP = 1.5

  type vector_i
     private
     integer          :: end  = 0
     integer          :: size = 0
     integer, pointer :: d(:) => NULL()
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

  subroutine vector_build_v_i(vi, values)
    type(vector_i), intent(out) :: vi
    integer,        intent(in ) :: values(:)
    if ( .not. associated(vi%d) ) then
       vi%size = size(values)
       allocate(vi%d(vi%size))
       vi%d     = values
       vi%end   = vi%size
    end if
  end subroutine vector_build_v_i

  subroutine vector_build_i(vi, size, value)
    type(vector_i), intent(out)           :: vi
    integer,        intent(in ), optional :: size
    integer,        intent(in ), optional :: value
    integer                               :: v
    v = 0
    if ( present(value) ) v      = value
    if ( present(size ) ) vi%size = size
    if ( .not. associated(vi%d)  ) then
       allocate(vi%d(vi%size))
       vi%d = v
    end if
  end subroutine vector_build_i

  subroutine vector_destroy_i(vi)
    type(vector_i), intent(inout) :: vi
    call clear(vi)
  end subroutine vector_destroy_i

  pure function vector_empty_i(vi) result(r)
    logical                    :: r
    type(vector_i), intent(in) :: vi
    r = vi%size == 0
  end function vector_empty_i

  pure function vector_size_i(vi) result(r)
    integer                    :: r
    type(vector_i), intent(in) :: vi
    r = vi%size
  end function vector_size_i

  pure function vector_capacity_i(vi) result(r)
    integer                    :: r
    type(vector_i), intent(in) :: vi
    if ( associated(vi%d)  ) then
       r = size(vi%d)
    else
       r = 0
    endif
  end function vector_capacity_i

  pure function vector_at_i(vi, i) result(r)
    integer                    :: r
    integer,        intent(in) :: i
    type(vector_i), intent(in) :: vi
    r = vi%d(i)
  end function vector_at_i

  subroutine vector_at_set_i(vi, i, v)
    integer                       :: v
    integer,        intent(in   ) :: i
    type(vector_i), intent(inout) :: vi
    vi%d(i) = v
  end subroutine vector_at_set_i

  function vector_dataptr_i(vi,lo,hi) result(r)
    integer, pointer                     :: r(:)
    type(vector_i), intent(in)           :: vi
    integer,        intent(in), optional :: lo
    integer,        intent(in), optional :: hi

    if ( present(lo) .AND. present(hi) ) then
       r => vi%d(lo:hi)
    else if ( present(lo) ) then
       r => vi%d(lo:vi%end)
    else if ( present(hi) ) then
       r => vi%d(1:hi)
    else
       r => vi%d(1:vi%end)
    end if

  end function vector_dataptr_i

  subroutine vector_resize_i(vi, size, value)
    type(vector_i), intent(inout) :: vi
    integer, intent(in) :: size
    integer, intent(in), optional ::  value
    integer, pointer :: np(:)
    integer :: v

    v = 0
    if ( present(value) ) v = value
    if ( size <= vector_capacity_i(vi) ) then
       vi%d(1:vi%size) = vi%d(1:vi%end)
       vi%d(vi%size+1:size) = v
    else
       allocate(np(size))
       np(1:vi%size) = vi%d(1:vi%end)
       np(vi%size+1:size) = v
       if ( associated(vi%d) ) deallocate(vi%d)
       vi%d => np
    end if
    vi%size  = size
    vi%end   = size

  end subroutine vector_resize_i

  subroutine vector_reserve_i(vi, size)
    type(vector_i), intent(inout) :: vi
    integer,        intent(in)    :: size
    integer, pointer              :: np(:)

    if ( size <= vector_capacity_i(vi) ) return

    allocate(np(size))
    np(1:vi%size) = vi%d(1:vi%size)
    if ( associated(vi%d) ) deallocate(vi%d)
    vi%d => np

  end subroutine vector_reserve_i

  subroutine vector_shrink_wrap_i(vi)
    type(vector_i), intent(inout) :: vi
    integer, pointer              :: np(:)

    if ( size(vi%d) == vi%size ) return

    allocate(np(vi%size))
    np(1:vi%size) = vi%d(1:vi%size)
    if ( associated(vi%d) ) deallocate(vi%d)
    vi%d => np
    vi%end   = vi%size

  end subroutine vector_shrink_wrap_i

  subroutine vector_push_back_i(vi,v)
    type(vector_i), intent(inout) :: vi
    integer,        intent(in   ) :: v

    if ( vi%size >= vector_capacity_i(vi) ) then
       call vector_reserve_i(vi,max(vi%size+4,int(BUMP*vi%size)))
    end if
    vi%size     = vi%size + 1
    vi%end      = vi%end  + 1
    vi%d(vi%end) = v

  end subroutine vector_push_back_i

  subroutine vector_pop_back_i(vi)
    type(vector_i), intent(inout) :: vi
    vi%size = vi%size - 1
    vi%end  = vi%end  - 1
  end subroutine vector_pop_back_i

  pure function vector_front_i(vi) result(r)
    type(vector_i), intent(in) :: vi
    integer :: r
    r = vi%d(1)
  end function vector_front_i

  pure function vector_back_i(vi) result(r)
    type(vector_i), intent(in) :: vi
    integer :: r
    r = vi%d(vi%end)
  end function vector_back_i

  subroutine vector_reverse_i(vi)
    type(vector_i), intent(inout) :: vi
    integer :: i,n

    n = size(vi)
    do i = 1, n/2
       call sw(vi%d(i),vi%d(n-i+1))
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

  function vector_erase_i(vi,i) result(r)
    type(vector_i), intent(inout) :: vi
    integer, intent(in) :: i
    integer :: r
    vi%d(i:vi%size-1) = vi%d(i+1:vi%size)
    vi%size = vi%size - 1
    vi%end  = vi%end  - 1
    r = i
  end function vector_erase_i

  subroutine vector_insert_i(vi, i, elem)
    type(vector_i), intent(inout) :: vi
    integer,        intent(in   ) :: i
    integer,        intent(in   ) :: elem

    if ( vi%size >= vector_capacity_i(vi) ) then
       call vector_reserve_i(vi,max(vi%size+4,int(BUMP*vi%size)))
    end if
    vi%size = vi%size + 1
    vi%end  = vi%end  + 1
    vi%d(i+1:vi%size) = vi%d(i:vi%size-1)
    vi%d(i) = elem

  end subroutine vector_insert_i

  subroutine vector_insert_n_i(vi, i, elem, n)
    type(vector_i), intent(inout) :: vi
    integer,        intent(in   ) :: n
    integer,        intent(in   ) :: i
    integer,        intent(in   ) :: elem
    integer j
    do j = 1, n
       call vector_insert_i(vi,i,elem)
    end do
  end subroutine vector_insert_n_i

  subroutine vector_clear_i(vi)
    type(vector_i), intent(inout) :: vi
    if ( associated(vi%d) ) then
       deallocate(vi%d)
       vi%d    => Null()
       vi%end  = 0
       vi%size = 0
    end if
  end subroutine vector_clear_i

  pure function vector_equal_i(vi1, vi2) result(r)
    logical                    :: r
    type(vector_i), intent(in) :: vi1, vi2
    integer                    :: i

    r = .TRUE.
    if ( size(vi1) .ne. size(vi2) ) then
       r = .FALSE.
       return
    end if
    do i = 1, size(vi1)
       if ( vi1%d(i) .NE. vi2%d(i) ) then
          r = .FALSE.
          exit
       end if
    end do
  end function vector_equal_i

  pure function vector_not_equal_i(vi1, vi2) result(r)
    logical :: r
    type(vector_i), intent(in) :: vi1, vi2
    r = .not. equal(vi1,vi2)
  end function vector_not_equal_i

   subroutine vector_print_i(vi, str, unit)
     use bl_IO_module
     type(vector_i), intent(in) :: vi
     character (len=*), intent(in), optional :: str
     integer, intent(in), optional :: unit
     integer :: un

     un = unit_stdout(unit)
     if ( present(str) ) write(unit=un, fmt='(A)', advance='no') str
     if ( empty(vi) ) then
        write(unit=un, fmt='("Empty")')
     else
        write(unit=un, fmt=*) vi%d(1:vi%end)
     end if

   end subroutine vector_print_i

end module vector_i_module
