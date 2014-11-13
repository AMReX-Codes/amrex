module sort_i_module
  ! Adapted from Meissner, Adams, et.al., and NR.

  use bl_types

  implicit none

  interface sort
     module procedure heapsort_i
     module procedure heapsort_indirect_i
  end interface sort

  interface stable_sort
     module procedure mergesort_i
     module procedure mergesort_indirect_i
     module procedure mergesort_indirect_ll
  end interface stable_sort

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

  subroutine mergesort_i( array, cmp )
    integer, dimension(:), intent(inout) :: array
    interface
       function cmp(a,b) result(r)         
         implicit none
         logical :: r
         integer, intent(in) :: a, b
       end function cmp
    end interface
    optional cmp

    integer, allocatable, dimension(:) :: tmp, tmp1, tmp2, tmp3, tmp4, tmpa, tmpb
    integer :: narray, lo(4), hi(4), nt(4), na, nb

    narray = size(array)

    if (present(cmp)) then
       if (narray .ge. 8) then

          lo(1) = 1
          lo(2) = 1+narray/4
          lo(3) = 1+narray/2
          lo(4) = 1+(narray*3)/4
          hi(1) = lo(2)-1
          hi(2) = lo(3)-1
          hi(3) = lo(4)-1
          hi(4) = narray
          nt = (hi-lo+2)/2
          allocate(tmp1(nt(1)))
          allocate(tmp2(nt(2)))
          allocate(tmp3(nt(3)))
          allocate(tmp4(nt(4)))
          na = hi(2)
          nb = narray-na
          allocate(tmpa(na))
          allocate(tmpb(nb))

          ! wz: I know this is ugly.  But I don't want to use nested parallel regions.
          !     I tried OMP TASK, but I couldn't get it work with Intel compiler.

          !$omp parallel 

          !$omp sections
          !$omp section
          call mg_sort_cmp(array(lo(1):hi(1)),tmp1)
          !$omp section
          call mg_sort_cmp(array(lo(2):hi(2)),tmp2)
          !$omp section
          call mg_sort_cmp(array(lo(3):hi(3)),tmp3)
          !$omp section
          call mg_sort_cmp(array(lo(4):hi(4)),tmp4)
          !$omp end sections

          !$omp sections
          !$omp section
          call merge2_cmp(array(lo(1):hi(1)),array(lo(2):hi(2)),tmpa)
          !$omp section
          call merge2_cmp(array(lo(3):hi(3)),array(lo(4):hi(4)),tmpb)
          !$omp end sections          

          !$omp end parallel

          call merge2_cmp(tmpa, tmpb, array)

       else
          allocate(tmp((narray+1)/2))
          call mg_sort_cmp(array, tmp)
       end if
    else
       allocate(tmp((narray+1)/2))
       call mg_sort_no_cmp(array, tmp)
    end if
  contains

    recursive subroutine mg_sort_cmp(a, tmp)
      integer, intent(inout) :: a(:) 
      integer                :: tmp(:)
      
      integer :: n, nh

      n = size(a)
      
      if (n .eq. 1) then

         return

      else if (n .eq. 2) then

         if (cmp(a(2),a(1))) then
            a(1:2) = a(2:1:-1)
         end if

      else

         nh = (n+1)/2

         call mg_sort_cmp(a(1:nh),tmp)
         call mg_sort_cmp(a(nh+1:),tmp)

         if (cmp(a(nh+1),a(nh))) then
            tmp(1:nh) = a(1:nh)
            call merge_cmp(tmp(1:nh), a)
         end if

      end if

      return
    end subroutine mg_sort_cmp

    subroutine merge_cmp(tmp, a)
      integer, intent(inout) :: a(:) 
      integer, intent(in   ) :: tmp(:) 
      
      integer :: n, nh, i, j, k

      n = size(a)
      nh = size(tmp)

      i=1 
      j=nh+1 
      k=1

      do while (i .le. nh .and. j .le. n)
         if (cmp(a(j),tmp(i))) then
            a(k) = a(j)
            j = j+1
         else
            a(k) = tmp(i)
            i = i+1
         end if
         k = k+1
      end do

      if (i .le. nh) a(k:) = tmp(i:nh)

      return
    end subroutine merge_cmp

    recursive subroutine mg_sort_no_cmp(a, tmp)
      integer, intent(inout) :: a(:) 
      integer                :: tmp(:)
      
      integer :: n, nh

      n = size(a)
      
      if (n .eq. 1) then

         return

      else if (n .eq. 2) then

         if (a(2) < a(1)) then
            a(1:2) = a(2:1:-1)
         end if

      else

         nh = (n+1)/2
         call mg_sort_no_cmp(a(1:nh),tmp)
         call mg_sort_no_cmp(a(nh+1:),tmp)
         
         if (a(nh+1) < a(nh)) then
            tmp(1:nh) = a(1:nh)
            call merge_no_cmp(tmp(1:nh), a)
         end if

      end if

      return
    end subroutine mg_sort_no_cmp

    subroutine merge_no_cmp(tmp, a)
      integer, intent(inout) :: a(:)
      integer, intent(in   ) :: tmp(:) 
      
      integer :: n, nh, i, j, k

      n = size(a)
      nh = size(tmp)

      i=1 
      j=nh+1 
      k=1

      do while (i .le. nh .and. j .le. n)
         if (a(j) < tmp(i)) then
            a(k) = a(j)
            j = j+1
         else
            a(k) = tmp(i)
            i = i+1
         end if
         k = k+1
      end do

      if (i .le. nh) a(k:) = tmp(i:nh)

      return
    end subroutine merge_no_cmp

    subroutine merge2_cmp(t1, t2, a)
      integer, intent(in   ) :: t1(:), t2(:)
      integer, intent(inout) :: a(:)
      integer :: i, j, k, n1, n2
      n1 = size(t1)
      n2 = size(t2)

      i = 1; j = 1; k = 1
      
      do while (i .le. n1 .and. j .le. n2)
         if (cmp(t2(j),t1(i))) then
            a(k) = t2(j)
            j = j+1
         else
            a(k) = t1(i)
            i = i+1
         end if
         k = k+1
      end do

      if (i .le. n1) then
         a(k:) = t1(i:)
      else
         a(k:) = t2(j:)
      end if
    end subroutine merge2_cmp

  end subroutine mergesort_i

  subroutine mergesort_indirect_i( array, iarray)
    integer, dimension(:), intent(in) :: array
    integer, dimension(:), intent(out) :: iarray

    integer :: tmp((size(array)+1)/2)
    integer :: i, n

    n = size (array)
    iarray = (/(i,i=1,n)/)
    
    call mg_sort_no_cmp(iarray, tmp)
    
  contains

    recursive subroutine mg_sort_no_cmp(idx, tmp)
      integer, intent(inout) :: idx(:) 
      integer                :: tmp(:)
      
      integer :: n, nh

      n = size(idx)
      
      if (n .eq. 1) then

         return

      else if (n .eq. 2) then

         if (array(idx(2)) < array(idx(1))) then
            idx(1:2) = idx(2:1:-1)
         end if

      else

         nh = (n+1)/2
         call mg_sort_no_cmp(idx(1:nh),tmp)
         call mg_sort_no_cmp(idx(nh+1:),tmp)
         
         if (array(idx(nh+1)) < array(idx(nh))) then
            tmp(1:nh) = idx(1:nh)
            call merge_no_cmp(tmp(1:nh), idx)
         end if

      end if

      return
    end subroutine mg_sort_no_cmp

    subroutine merge_no_cmp(tmp, idx)
      integer, intent(inout) :: idx(:)
      integer, intent(in   ) :: tmp(:) 
      
      integer :: n, nh, i, j, k

      n = size(idx)
      nh = size(tmp)

      i=1 
      j=nh+1 
      k=1

      do while (i .le. nh .and. j .le. n)
         if (array(idx(j)) < array(tmp(i))) then
            idx(k) = idx(j)
            j = j+1
         else
            idx(k) = tmp(i)
            i = i+1
         end if
         k = k+1
      end do

      if (i .le. nh) idx(k:) = tmp(i:nh)

      return
    end subroutine merge_no_cmp

  end subroutine mergesort_indirect_i

  subroutine mergesort_indirect_ll( array, iarray, ascending)
    integer(kind=ll_t), dimension(:), intent(in) :: array
    integer           , dimension(:), intent(out) :: iarray
    logical, intent(in), optional :: ascending

    logical :: lascending
    integer :: tmp((size(array)+1)/2)
    integer :: i, n

    lascending = .true. ;  if (present(ascending)) lascending = ascending

    n = size (array)
    iarray = (/(i,i=1,n)/)
    
    call mg_sort_no_cmp(iarray, tmp)

  contains

    recursive subroutine mg_sort_no_cmp(idx, tmp)
      integer, intent(inout) :: idx(:) 
      integer                :: tmp(:)
      
      integer :: n, nh

      n = size(idx)
      
      if (n .eq. 1) then
         
         return
         
      else if (n .eq. 2) then

         if (      lascending .and. array(idx(2)) < array(idx(1)) .or. &
              .not.lascending .and. array(idx(2)) > array(idx(1)) ) then
            idx(1:2) = idx(2:1:-1)
         end if

      else

         nh = (n+1)/2
         call mg_sort_no_cmp(idx(1:nh),tmp)
         call mg_sort_no_cmp(idx(nh+1:),tmp)

         if (      lascending .and. array(idx(nh+1)) < array(idx(nh)) .or. &
              .not.lascending .and. array(idx(nh+1)) > array(idx(nh)) ) then
            tmp(1:nh) = idx(1:nh)
            call merge_no_cmp(tmp(1:nh), idx)
         end if

      end if

      return
    end subroutine mg_sort_no_cmp

    subroutine merge_no_cmp(tmp, idx)
      integer, intent(inout) :: idx(:)
      integer, intent(in   ) :: tmp(:) 
      
      integer :: n, nh, i, j, k

      n = size(idx)
      nh = size(tmp)

      i=1 
      j=nh+1 
      k=1

      do while (i .le. nh .and. j .le. n)
         if (      lascending .and. array(idx(j)) < array(tmp(i)) .or. &
              .not.lascending .and. array(idx(j)) > array(tmp(i)) ) then
            idx(k) = idx(j)
            j = j+1
         else
            idx(k) = tmp(i)
            i = i+1
         end if
         k = k+1
      end do
      
      if (i .le. nh) idx(k:) = tmp(i:nh)

      return
    end subroutine merge_no_cmp

  end subroutine mergesort_indirect_ll

end module sort_i_module
