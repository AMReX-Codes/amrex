!! Uniform method of printing error/warning messages.  Error messages
!! print a message to the default output unit and then call the PARALLEL_ABORT
!! to the process (or set of processes if parallel)

module bl_error_module

  use bl_types

  implicit none

  !! Print an error message consisting of text followed by an optional
  !! scalar variable of type character(len=*), integer, or real; the
  !! process terminates.
  interface bl_error
     module procedure bl_error0
     module procedure bl_error1_ch
     module procedure bl_error1_i
     module procedure bl_error1_d
     module procedure bl_error1_s
     module procedure bl_error1_z
     module procedure bl_error1_c
  end interface

  !! Print an warning message consisting of text followed by an optional
  !! scalar variable of type character(len=*), integer, or real.
  interface bl_warn
     module procedure bl_warn0
     module procedure bl_warn1_ch
     module procedure bl_warn1_i
     module procedure bl_warn1_d
     module procedure bl_warn1_s
     module procedure bl_warn1_z
     module procedure bl_warn1_c
  end interface

  !! Print an error message consisting of text when logical condition
  !! in arguments are not all .TRUE.
  interface bl_assert
     module procedure bl_assert1
     module procedure bl_assert2
     module procedure bl_assert3
     module procedure bl_assert4
     module procedure bl_assert_v
  end interface

contains

  subroutine bl_error0(str)
    use parallel
    character(len=*), intent(in), optional :: str
    if ( present(str) ) then
       write(*,*) "BOXLIB ERROR: ", str
    else
       write(*,*) "BOXLIB ERROR"
    end if
    call parallel_abort()
  end subroutine bl_error0

  subroutine bl_error1_ch(str, str1)
    use parallel
    character(len=*), intent(in) :: str, str1
    write(*,fmt=*) "BOXLIB ERROR: ", str, str1
    call parallel_abort()
  end subroutine bl_error1_ch

  subroutine bl_error1_i(str, val)
    use parallel
    character(len=*), intent(in) :: str
    integer, intent(in) :: val
    write(*,fmt=*) "BOXLIB ERROR: ", str, val
    call parallel_abort()
  end subroutine bl_error1_i

  subroutine bl_error1_d(str, val)
    use parallel
    character(len=*), intent(in) :: str
    real(kind=dp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB ERROR: ", str, val
    call parallel_abort()
  end subroutine bl_error1_d

  subroutine bl_error1_s(str, val)
    use parallel
    character(len=*), intent(in) :: str
    real(kind=sp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB ERROR: ", str, val
    call parallel_abort()
  end subroutine bl_error1_s

  subroutine bl_error1_z(str, val)
    use parallel
    character(len=*), intent(in) :: str
    complex(kind=dp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB ERROR: ", str, val
    call parallel_abort()
  end subroutine bl_error1_z

  subroutine bl_error1_c(str, val)
    use parallel
    character(len=*), intent(in) :: str
    complex(kind=sp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB ERROR: ", str, val
    call parallel_abort()
  end subroutine bl_error1_c

  subroutine bl_warn0(str)
    character(len=*), intent(in) :: str
    write(*,fmt=*) "BOXLIB WARN: ", str
  end subroutine bl_warn0

  subroutine bl_warn1_ch(str, str1)
    character(len=*), intent(in) :: str, str1
    write(*,fmt=*) "BOXLIB WARN: ", str, str1
  end subroutine bl_warn1_ch

  subroutine bl_warn1_i(str, val)
    character(len=*), intent(in) :: str
    integer, intent(in) :: val
    write(*,fmt=*) "BOXLIB WARN: ", str, val
  end subroutine bl_warn1_i

  subroutine bl_warn1_d(str, val)
    character(len=*), intent(in) :: str
    real(kind=dp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB WARN: ", str, val
  end subroutine bl_warn1_d

  subroutine bl_warn1_s(str, val)
    character(len=*), intent(in) :: str
    real(kind=sp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB WARN: ", str, val
  end subroutine bl_warn1_s

  subroutine bl_warn1_z(str, val)
    character(len=*), intent(in) :: str
    complex(kind=dp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB WARN: ", str, val
  end subroutine bl_warn1_z

  subroutine bl_warn1_c(str, val)
    character(len=*), intent(in) :: str
    complex(kind=sp_t), intent(in) :: val
    write(*,fmt=*) "BOXLIB WARN: ", str, val
  end subroutine bl_warn1_c

  !! Stolen from Numerical Recepies
  !! If COND is true, nothing; if COND is false, call BL_ERROR_CH, which
  !! terminates the process

  subroutine bl_assert1(n1, str)
    character(len=*), intent(in) :: str
    logical, intent(in) :: n1
    if ( .not. n1 ) then
       call bl_error("ASSERTION FAILED:1: ", str)
    end if
  end subroutine bl_assert1

  subroutine bl_assert2(n1, n2, str)
    character(len=*), intent(in) :: str
    logical, intent(in) :: n1, n2
    if ( .not. (n1 .and. n2) ) then
       call bl_error("ASSERTION FAILED:2: ", str)
    end if
  end subroutine bl_assert2

  subroutine bl_assert3(n1, n2, n3, str)
    character(len=*), intent(in) :: str
    logical, intent(in) :: n1, n2, n3
    if ( .not. (n1 .and. n2 .and. n3) ) then
       call bl_error("ASSERTION FAILED:3: ", str)
    end if
  end subroutine bl_assert3

  subroutine bl_assert4(n1, n2, n3, n4, str)
    character(len=*), intent(in) :: str
    logical, intent(in) :: n1, n2, n3, n4
    if ( .not. (n1 .and. n2 .and. n3 .and. n4) ) then
       call bl_error("ASSERTION FAILED:4: ", str)
    end if
  end subroutine bl_assert4

  subroutine bl_assert_v(n, str)
    character(len=*), intent(in) :: str
    logical, dimension(:), intent(in) :: n
    if ( .not. all(n) ) then
       call bl_error("ASSERTION FAILED:v: ", str)
    end if
  end subroutine bl_assert_v

end module bl_error_module
