module list_box_module
  use box_module;  
  implicit none

  type list_box_node
     private
     type(box) :: v ! = 
     type(list_box_node), pointer :: prev => NULL()
     type(list_box_node), pointer :: next => NULL()
  end type list_box_node

  type list_box
     private
     integer :: size = 0
     type(list_box_node), pointer:: head => NULL()
     type(list_box_node), pointer:: tail => NULL()
  end type list_box

  interface value
     module procedure list_node_value_box
  end interface value
  interface get
     module procedure list_node_get_box
  end interface get
  interface set
     module procedure list_node_set_box
  end interface set

  interface assignment (=)
     module procedure list_node_set_box
     module procedure list_node_get_box
  end interface

  interface next
     module procedure list_node_next_box
  end interface next

  interface prev
     module procedure list_node_prev_box
  end interface prev

  interface build
     module procedure list_build_v_box
     module procedure list_build_sn_box
     module procedure list_build_box
  end interface build

  interface destroy
     module procedure list_destroy_box
  end interface destroy

  interface size
     module procedure list_size_box
  end interface size

  interface empty
     module procedure list_empty_box
  end interface empty

  interface begin
     module procedure list_begin_box
  end interface begin

  interface end
     module procedure list_end_box
  end interface end

  interface back
     module procedure list_back_box
  end interface back

  interface front
     module procedure list_front_box
  end interface front

  interface push_back
     module procedure list_push_back_box
  end interface push_back

  interface insert
     module procedure list_insert_box
  end interface insert

  interface push_front
     module procedure list_push_front_box
  end interface push_front

  interface pop_back
     module procedure list_pop_back_box
  end interface pop_back

  interface pop_front
     module procedure list_pop_front_box
  end interface pop_front

  interface clear
     module procedure list_clear_box
  end interface clear

  interface erase
     module procedure list_erase_box
     module procedure list_erase_range_box
  end interface erase

  interface splice
     module procedure list_splice_box
     ! module procedure list_splice_v_box
  end interface splice

  interface unique
     module procedure list_unique_box
  end interface unique

  interface remove
     module procedure list_remove_box
  end interface remove

  interface remove_if
     module procedure list_remove_if_box
  end interface remove_if

  interface reverse
     module procedure list_reverse_box
  end interface reverse

  interface swap
     module procedure list_swap_box
  end interface swap

  interface operator (.EQ.)
     module procedure list_equal_box
  end interface

  interface operator (.NE.)
     module procedure list_equal_box
  end interface

contains


  ! Node accessors.

  function list_node_value_box(n) result(r)
    type(list_box_node), intent(in) :: n
    type(box) :: r

    r = n%v

  end function list_node_value_box

  subroutine list_node_set_box(n, v)
    type(list_box_node), intent(inout) :: n
    type(box), intent(in) ::  v

    n%v = v

  end subroutine list_node_set_box

  subroutine list_node_get_box(v, n)
    type(list_box_node), intent(in) :: n
    type(box), intent(out) ::  v

    v = n%v

  end subroutine list_node_get_box

  function list_node_next_box(n) result(r)
    type(list_box_node), pointer :: n, r

    r => n%next

  end function list_node_next_box

  function list_node_prev_box(n) result(r)
    type(list_box_node), pointer :: n, r

    r => n%prev

  end function list_node_prev_box

  subroutine list_build_box(r)
    type(list_box), intent(out) :: r

    r%size = 0

  end subroutine list_build_box

  subroutine list_build_v_box(r, d)
    type(list_box), intent(out) :: r
    type(box), dimension(:), intent(in) :: d
    integer :: i

    do i = 1, size(d)
       call push_back(r, d(i))
    end do
    r%size = size(d)

  end subroutine list_build_v_box

  subroutine list_build_sn_box(r, n, v)
    type(list_box), intent(out) :: r
    integer, intent(in) :: n
    type(box), intent(in) :: v
    integer i

    do i = 1, n
       call push_back(r, v)
    end do
    r%size = n

  end subroutine list_build_sn_box

  subroutine list_destroy_box(r)
    type(list_box), intent(inout) :: r

    call clear(r)
    r%size = 0

  end subroutine list_destroy_box

  subroutine list_copy_box(l1,l2)
    type(list_box), intent(inout) :: l1
    type(list_box), intent(in) :: l2
    type(list_box_node), pointer :: b

    call clear(l1)
    b => begin(l2)
    do while ( associated(b) )
       call push_back(l1, value(b))
       b => next(b)
    end do
    l1%size = l2%size

  end subroutine list_copy_box

  function list_size_box(l) result (r)
    type(list_box), intent(in) :: l
    integer :: r

    r = l%size

  end function list_size_box

  function list_begin_box(l) result(r)
    type(list_box), intent(in) :: l
    type(list_box_node), pointer :: r

    r => l%head

  end function list_begin_box

  function list_end_box(l) result(r)
    type(list_box), intent(in) :: l
    type(list_box_node), pointer :: r

    r => l%tail%next

  end function list_end_box

  function list_empty_box(l) result (r)
    type(list_box), intent(in) :: l
    logical:: r

    r = .not.(associated(l%head) .AND. associated(l%tail))

  end function list_empty_box

  function list_find_box(l, number) result ( r )
    type(list_box), intent(in) :: l
    type(box), intent(in) :: number
    type(list_box_node), pointer :: r

    r => l%head
    do
       if ( .not. associated(r) ) return
       if ( number == value(r) ) return
       r => r%next
    end do

  end function list_find_box

  function list_equal_box(l1,l2) result (r)
    type(list_box), intent(in) :: l1, l2
    logical :: r
    type(list_box_node), pointer :: p1, p2

    p1 => l1%head
    p2 => l2%head
    do while( associated(p1) .AND. associated(p2) )
       if ( value(p1) /= value(p2) ) then
          r = .FALSE.
          return
       end if
       p1 => p1%next
       p2 => P2%next
    end do
    if ( .not. associated(p1) .AND. .not. associated(p2) ) then
       r = .TRUE.
    else
       r = .FALSE.
    end if

  end function list_equal_box

  function list_not_equal_box(l1,l2) result (r)
    type(list_box), intent(in) :: l1, l2
    logical :: r

    r = .not. list_equal_box(l1,l2)

  end function list_not_equal_box

  subroutine list_push_front_box(l,v)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(box), intent(in) :: v
    type(list_box_node), pointer :: n
    integer :: ierr
    allocate(n, stat = ierr)
    if ( ierr /= 0 ) then
       call bl_error("list_push_front_box: failed to allocate memory")
    end if
    n%v = v
    n%next => l%head
    l%head => n
    if ( .not. associated(l%tail) ) then
       l%tail => l%head
    else
       l%head%next%prev => n
    end if
    l%size = l%size + 1

  end subroutine list_push_front_box

  subroutine list_push_back_box(l,v)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(box), intent(in) :: v
    type(list_box_node), pointer :: n
    integer :: ierr
    allocate(n, stat = ierr)
    if ( ierr /= 0 ) then
       call bl_error("list_push_back_box: failed to allocate memory")
    end if
    n%v = v
    n%prev => l%tail
    l%tail => n
    if ( .not. associated (l%head ) ) then
       l%head => l%tail
    else
       l%tail%prev%next => n
    end if
    l%size = l%size + 1

  end subroutine list_push_back_box

  subroutine list_insert_box(l, ln, v)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: ln, n
    type(box), intent(in) :: v
    integer :: ierr

    if ( associated(ln, l%head) ) then
       call push_front(l, v)
    else
       allocate(n, stat = ierr)
       if ( ierr /= 0 ) then
          call bl_error("list_push_back_box: failed to allocate memory")
       end if
       n%v = v
       n%prev => ln%prev
       n%next => ln
       ln%prev%next => n
       ln%prev => n
    end if
    l%size = l%size + 1

  end subroutine list_insert_box

  function list_front_box(l) result ( r )
    type(list_box), intent(in) :: l
    type(box) :: r

    if ( empty(l) ) STOP 'FRONT: on empty list'
    r = value(l%head)

  end function list_front_box

  function list_back_box(l) result ( r )
    type(list_box), intent(in) :: l
    type(box) :: r

    if ( empty(l) ) STOP 'BACK: on empty list'
    r = value(l%tail)

  end function list_back_box

  subroutine list_pop_front_box(l)
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: h, p

    if ( empty(l) ) STOP 'POP_FRONT: on empty list'
    h => l%head
    p => erase(l, h)

  end subroutine list_pop_front_box

  subroutine list_pop_back_box(l)
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: b, p

    if ( empty(l) ) STOP 'POP_BACK: on empty list'
    b => l%tail
    p => erase(l, b)

  end subroutine list_pop_back_box

  function list_erase_box(l, beg) result (r)
    use bl_error_module
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: beg, r
    integer :: ierr

    if ( associated(beg, l%tail) ) then
       r => Null()
    else
       r => beg%next
    end if
    if ( associated(beg, l%head) .AND. associated(beg, l%tail) ) then
       l%head => Null()
       l%tail => Null()
    else if ( associated(beg, l%head) ) then
       l%head => beg%next
       l%head%prev => Null()
    else if ( associated(beg, l%tail) ) then
       l%tail => beg%prev
       l%tail%next => Null()
    else
       beg%next%prev => beg%prev
       beg%prev%next => beg%next
    end if
    deallocate(beg, stat = ierr)
    if ( ierr /= 0 ) then
       call bl_error("list_erase_box: failed to deallocate")
    end if
    l%size = l%size - 1

  end function list_erase_box

  function list_erase_range_box(l, beg, end) result(rr)
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: beg, end
    type(list_box_node), pointer :: r, e, rr

    r => beg
    e => end
    if ( .not. associated(e) ) then
       do while ( associated(r) )
          r => erase(l,r)
       end do
    else
       do while ( .not. associated(r, e) )
          r => erase(l,r)
       end do
    end if
    rr => r

  end function list_erase_range_box

  subroutine list_clear_box(l)
    type(list_box), intent(inout) :: l
    type(list_box_node), pointer :: r

    r => erase(l, l%head, Null(r))
    if ( .not. empty(l) ) STOP 'CLEAR: is broken'
    if ( l%size /= 0 ) stop 'CLEAR: size is broken'

  end subroutine list_clear_box

  subroutine list_unique_box(l, tst)
    interface
       function tst(p,q) result(r)
         use box_module;  
         implicit none
         logical :: r
         type(box), intent(in) :: p, q
       end function tst
    end interface
    optional tst
    type(list_box), intent(inout) :: l

    if ( present(tst) ) then
       call un(tst)
    else
       call un(box_equal)
    end if

  contains

    subroutine un(tst)
      use bl_error_module
      interface
         function tst(p,q) result(r)
           use box_module;  
           implicit none
           logical :: r
           type(box), intent(in) :: p, q
         end function tst
      end interface
      type(list_box_node), pointer :: head, tail, e, p
      integer :: ierr

      if ( size(l) <= 1 ) return
      l%size = 1
      head => l%head
      tail => head
      e => head%next
      do while ( associated(e) )
         if ( tst(tail%v, e%v) ) then
            p => e
            e => e%next
            deallocate(p, stat=ierr)
            if ( Ierr /= 0 ) then
               call bl_error("list_unique_op_box: failed to deallocated memory")
            end if
            cycle
         end if
         tail%next => e
         e%prev = tail
         tail => e
         e => tail%next
         l%size = l%size + 1
      end do
      l%head => head
      l%head%prev => Null()
      l%tail => tail
      l%tail%next => Null()

    end subroutine un

  end subroutine list_unique_box

  subroutine list_splice_box(l1, l2)
    type(list_box), intent(inout) :: l1, l2

    if ( empty(l2) ) return
    if ( empty(l1) ) then
       l1%head => l2%head
       l1%tail => l2%tail
    else
       l1%tail%next => l2%head
       l2%head%prev => l1%tail
       l1%tail => l2%tail
    end if
    nullify(l2%head, l2%tail)
    l1%size = l2%size + l1%size
    l2%size = 0

  end subroutine list_splice_box

  subroutine list_splice_v_box(l1, pos, l2, first, last)
    type(list_box), intent(inout) :: l1, l2
    type(list_box_node), pointer :: pos
    type(list_box_node), pointer, optional :: first,last
    type(list_box_node), pointer :: p,f,l

    if ( empty(l2) ) return
    stop 'SPLICE_V: NOT WRITTEN YET'
    p => pos
    if ( present(first) ) then
       f => first
    else
       f => begin(l2)
    end if
    if ( present(last) ) then
       l => last
    else
       l => f
    end if
    if ( associated(l) ) then
       do while ( .not. associated(f,l) )
          call splice_one(p,f)
          f => next(f)
       end do
    else
       do while ( associated(f) )
          call splice_one(p, f)
          f => next(f)
       end do
    end if

  contains

    subroutine splice_one(p,f)
      type(list_box_node), pointer :: p, f

      ! remove f from its list
      f%prev%next => f%next
      f%next%prev => f%prev
      ! place f into list ahead of p
      f%next => p
      f%prev => p%prev
      ! place f ahead of p
      p%prev%next => f
      p%prev => f
      l1%size = l1%size + 1
      l2%size = l2%size - 1

    end subroutine splice_one

  end subroutine list_splice_v_box


  subroutine list_reverse_box(l)
    type(list_box), intent(inout) :: l
    type(list_box) :: l1

    call build(l1)
    do while ( associated(begin(l)) )
       call push_front(l1, front(l))
       call pop_front(l)
    end do
    call destroy(l)
    l = l1

  end subroutine list_reverse_box

  subroutine list_remove_box(l, val)
    type(list_box), intent(inout) :: l
    type(box), intent(in) :: val
    type(list_box_node), pointer :: p

    p => begin(l)
    do while ( associated(p) )
       if ( value(p) == val ) then
          p => erase(l,p)
       else
          p => next(p)
       end if
    end do

  end subroutine list_remove_box

  subroutine list_remove_if_box(l, tst)
    type(list_box), intent(inout) :: l
    interface
       function tst(p) result(r)
         use box_module;  
         implicit none
         logical :: r
         type(box), intent(in) :: p
       end function tst
    end interface
    type(list_box_node), pointer :: p

    p => begin(l)
    do while ( associated(p) )
       if ( tst(value(p)) ) then
          p => erase(l, p)
       else
          p => next(p)
       end if
    end do

  end subroutine list_remove_if_box

  subroutine list_swap_box(l1,l2)
    type(list_box), intent(inout) :: l1, l2
    type(list_box) :: lt

    lt = l1; l1 = l2; l2 = lt

  end subroutine list_swap_box


end module list_box_module
