module sort_box_module
  ! Adapted from Meissner, Adams, et.al., and NR.

  use box_module

  implicit none

  interface box_sort
     module procedure heapsort_box
     module procedure heapsort_indirect_box
  end interface box_sort

contains

  subroutine heapsort_box( array, cmp)
    type(box), dimension(:), intent(in out) :: array
    interface
       function cmp(a,b) result(r)
         use box_module
         implicit none
         logical :: r
         type(box), intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    if ( present(cmp) ) then
       call he_sort_cmp(array,cmp)
    else
       call he_sort_no_cmp(array)
    end if

  contains

    subroutine box_swap(a,b)
      type(box), intent(inout) :: a, b
      type(box) :: t
      t = a; a = b; b = t
    end subroutine box_swap

    subroutine he_sort_cmp(array,cmp)
      type(box), dimension(:), intent(in out) :: array
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_cmp(array,i,n,cmp)
      end do
      do i = n, 2, -1
         call box_swap(array(1),array(i))
         call peck_cmp(array,1, i-1,cmp)
      end do

    end subroutine he_sort_cmp

    subroutine he_sort_no_cmp(array)
      type(box), dimension(:), intent(in out) :: array
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_no_cmp(array,i,n)
      end do
      do i = n, 2, -1
         call box_swap(array(1),array(i))
         call peck_no_cmp(array,1, i-1)
      end do

    end subroutine he_sort_no_cmp

    subroutine peck_cmp(array,l,r,cmp)
      type(box), dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: i
      type(box) :: next

      next = array(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( cmp(array(i), array(i+1)) ) i = i+1
         end if
         if ( .not. cmp( next, array(i)) ) exit
         array(i/2) = array(i)
         i = 2*i
      end do
      array(i/2) = next

    end subroutine peck_cmp

    subroutine peck_no_cmp(array,l,r)
      type(box), dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      integer :: i

      type(box) :: next
      next = array(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( box_less(array(i), array(i+1)) ) i = i+1
         end if
         if ( .not. box_less( next, array(i)) ) exit
         array(i/2) = array(i)
         i = 2*i
      end do
      array(i/2) = next

    end subroutine peck_no_cmp

  end subroutine heapsort_box

  subroutine heapsort_indirect_box( array, iarray, cmp)
    type(box), dimension(:), intent(in) :: array
    integer, dimension(:), intent(out) :: iarray
    interface
       function cmp(a,b) result(r)
         use box_module
         implicit none
         logical :: r
         type(box), intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    if ( present(cmp) ) then
       call he_sort_cmp(iarray, cmp)
    else
       call he_sort_no_cmp(iarray)
    end if

  contains

    subroutine integer_swap(a,b)
      integer, intent(inout) :: a, b
      integer :: t
      t = a; a = b; b = t
    end subroutine integer_swap

    subroutine he_sort_no_cmp(iarray)
      integer, dimension(:), intent(out) :: iarray
      integer :: n, i

      n = size (array)
      iarray = (/(i,i=1,n)/)
      do i = n / 2, 1, -1
         call peck_no_cmp(iarray, i, n)
      end do
      do i = n, 2, -1
         call integer_swap(iarray(1), iarray(i))
         call peck_no_cmp(iarray, 1, i-1)
      end do

    end subroutine he_sort_no_cmp

    subroutine peck_no_cmp(iarray, l, r)
      integer, dimension(:), intent(inout) :: iarray
      integer, intent(in) :: r
      integer, intent(in) :: l
      integer :: i
      integer :: next

      next = iarray(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( box_less(array(iarray(i)), array(iarray(i+1))) ) i = i+1
         end if
         if ( .not. box_less( array(next), array(iarray(i))) ) exit
         iarray(i/2) = iarray(i)
         i = 2*i
      end do
      iarray(i/2) = next

    end subroutine peck_no_cmp

    subroutine he_sort_cmp(iarray, cmp)
      integer, dimension(:), intent(out) :: iarray
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: n, i

      n = size (array)
      iarray = (/(i,i=1,n)/)
      do i = n / 2, 1, -1
         call peck_cmp(iarray, i, n, cmp)
      end do
      do i = n, 2, -1
         call integer_swap(iarray(1), iarray(i))
         call peck_cmp(iarray, 1, i-1, cmp)
      end do

    end subroutine he_sort_cmp

    subroutine peck_cmp(iarray, l, r, cmp)
      integer, dimension(:), intent(inout) :: iarray
      integer, intent(in) :: r
      integer, intent(in) :: l
      interface
         function cmp(a,b) result(r)
           use box_module
           implicit none
           logical :: r
           type(box), intent(in) :: a, b
         end function cmp
      end interface
      integer :: i
      integer :: next

      next = iarray(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( cmp(array(iarray(i)), array(iarray(i+1))) ) i = i+1
         end if
         if ( .not. cmp( array(next), array(iarray(i))) ) exit
         iarray(i/2) = iarray(i)
         i = 2*i
      end do
      iarray(i/2) = next

    end subroutine peck_cmp

  end subroutine heapsort_indirect_box

end module sort_box_module

