module sort_i_module
  ! Adapted from Meissner, Adams, et.al., and NR.
  
  implicit none

  interface sort
     module procedure heapsort_i
     module procedure heapsort_indirect_i
  end interface sort

  private less_i
  private swap_i

contains

  pure function less_i(a,b) result(r)
    logical :: r
    integer, intent(in) :: a, b
    r = a < b
  end function less_i

  subroutine swap_i(a,b)
    integer, intent(inout) :: a, b
    integer :: t
    t = a; a = b; b = t
  end subroutine swap_i

  subroutine heapsort_i( array, cmp)
    integer, dimension(:), intent(in out) :: array
    interface
       function cmp(a,b) result(r)
         
         implicit none
         logical :: r
         integer, intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    if ( present(cmp) ) then
       call he_sort_cmp(array,cmp)
    else
       call he_sort_no_cmp(array)
    end if

  contains

    subroutine he_sort_cmp(array,cmp)
      integer, dimension(:), intent(in out) :: array
      interface
         function cmp(a,b) result(r)
           
           implicit none
           logical :: r
           integer, intent(in) :: a, b
         end function cmp
      end interface
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_cmp(array,i,n,cmp)
      end do
      do i = n, 2, -1
         call swap_i(array(1),array(i))
         call peck_cmp(array,1, i-1,cmp)
      end do

    end subroutine he_sort_cmp

    subroutine he_sort_no_cmp(array)
      integer, dimension(:), intent(in out) :: array
      integer :: n, i

      n = size (array)
      do i = n / 2, 1, -1
         call peck_no_cmp(array,i,n)
      end do
      do i = n, 2, -1
         call swap_i(array(1),array(i))
         call peck_no_cmp(array,1, i-1)
      end do

    end subroutine he_sort_no_cmp

    subroutine peck_cmp(array,l,r,cmp)
      integer, dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      interface
         function cmp(a,b) result(r)
           
           implicit none
           logical :: r
           integer, intent(in) :: a, b
         end function cmp
      end interface
      integer :: i
      integer :: next

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
      integer, dimension(:), intent(inout) :: array
      integer, intent(in) :: r
      integer, intent(in) :: l
      integer :: i

      integer :: next
      next = array(l)
      i = 2*l
      do while ( i <= r )
         if ( i < r ) then
            if ( less_i(array(i), array(i+1)) ) i = i+1
         end if
         if ( .not. less_i( next, array(i)) ) exit
         array(i/2) = array(i)
         i = 2*i
      end do
      array(i/2) = next

    end subroutine peck_no_cmp

  end subroutine heapsort_i

  subroutine heapsort_indirect_i( array, iarray, cmp)
    integer, dimension(:), intent(in) :: array
    integer, dimension(:), intent(out) :: iarray
    interface
       function cmp(a,b) result(r)
         implicit none
         logical :: r
         integer, intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    if ( present(cmp) ) then
       call he_sort_cmp(iarray, cmp)
    else
       call he_sort_no_cmp(iarray)
    end if

  contains

    subroutine he_sort_no_cmp(iarray)
      integer, dimension(:), intent(out) :: iarray
      integer :: n, i

      n = size (array)
      iarray = (/(i,i=1,n)/)
      do i = n / 2, 1, -1
         call peck_no_cmp(iarray, i, n)
      end do
      do i = n, 2, -1
         call swap_i(iarray(1), iarray(i))
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
            if ( less_i(array(iarray(i)), array(iarray(i+1))) ) i = i+1
         end if
         if ( .not. less_i( array(next), array(iarray(i))) ) exit
         iarray(i/2) = iarray(i)
         i = 2*i
      end do
      iarray(i/2) = next

    end subroutine peck_no_cmp

    subroutine he_sort_cmp(iarray, cmp)
      integer, dimension(:), intent(out) :: iarray
      interface
         function cmp(a,b) result(r)
           
           implicit none
           logical :: r
           integer, intent(in) :: a, b
         end function cmp
      end interface
      integer :: n, i

      n = size (array)
      iarray = (/(i,i=1,n)/)
      do i = n / 2, 1, -1
         call peck_cmp(iarray, i, n, cmp)
      end do
      do i = n, 2, -1
         call swap_i(iarray(1), iarray(i))
         call peck_cmp(iarray, 1, i-1, cmp)
      end do

    end subroutine he_sort_cmp

    subroutine peck_cmp(iarray, l, r, cmp)
      integer, dimension(:), intent(inout) :: iarray
      integer, intent(in) :: r
      integer, intent(in) :: l
      interface
         function cmp(a,b) result(r)
           implicit none
           logical :: r
           integer, intent(in) :: a, b
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

  end subroutine heapsort_indirect_i
end module sort_i_module
