
!=======================================================================
!
!     String handling routines:
!     Strings are handled differently in C++ and in FORTRAN.  In order
!     to simplify the framework strings are passed from FORTRAN to C++
!     as arrays of integer characters, terminated by the EOS symbol
!     which we set to -1
!     bl_str2int converts a FORTRAN string to an integer array,
!     bl_int2str converts an integer array to a FORTRAN string.
!      
!-----------------------------------------------------------------------

  subroutine bl_str2int(iarr, n, str)

    character*(*) :: str
    integer :: n, i, j
    integer :: iarr(n)
    integer, parameter :: EOS = -1

    if ( n .LE. len(str) ) then
       call bl_abort("bl_str2int: str to large for iarr")
    end if

    ! Make sure that IARR is empty
    do J = 1, N
       iarr(J) = ichar(' ')
    end do
    j = 1
    do i = 1, len(str)
       iarr(j) = ichar(str(i:i))
       j = j + 1
    end do

    ! EOS
    iarr(j) = EOS

  end subroutine bl_str2int

!-----------------------------------------------------------------------

  subroutine bl_int2str(str, iarr, n)

    character*(*) :: str
    integer :: n
    integer :: iarr(n)
    integer, parameter :: EOS = -1
    integer :: i

    do i = 1, LEN(str)
       str(i:i) = ' '
    end do
    do i = 1, n
       if ( i .GT. LEN(str) ) then
          call bl_abort("bl_int2str: iarr to large for str")
       end if
       if ( iarr(i) .EQ. EOS ) GO TO 100
       str(i:i) = char(iarr(i))
    end do

100 CONTINUE

  end subroutine bl_int2str
